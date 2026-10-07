"""Local sequence descriptors and strictly bound optional predictor evidence.

These descriptors are not SNIPR outcome predictors. No sequences leave the host.
"""
import hashlib
import math
from collections import Counter

from .common import digest, sequence

KD = dict(zip("ACDEFGHIKLMNPQRSTVWY", (1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8,
                                            1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3)))
PARAMETERS = {"hydropathy_window": 19, "hydropathy_alert": 1.6,
              "entropy_window": 12, "entropy_alert": 2.2,
              "stp_window": 30, "stp_alert": 0.65}


def sequence_hash(seq):
    return hashlib.sha256(seq.encode("ascii")).hexdigest()


def _segments(windows, width, minimum=False):
    result = []
    for start, value in windows:
        end = start + width - 1
        if result and start <= result[-1]["end"] + 1:
            result[-1]["end"] = max(end, result[-1]["end"])
            result[-1]["extreme_value"] = (min if minimum else max)(value, result[-1]["extreme_value"])
        else:
            result.append({"start": start, "end": end, "extreme_value": value})
    return result


def descriptors(seq, use_biopython=False):
    seq = sequence(seq)
    n = len(seq)
    kd = [sum(KD[a] for a in seq[i:i+19])/19 for i in range(max(0, n-18))]
    entropies = []
    for i in range(max(0, n-11)):
        counts = Counter(seq[i:i+12])
        entropies.append(-sum(v/12 * math.log2(v/12) for v in counts.values()))
    stp = [sum(a in "STP" for a in seq[i:i+30])/30 for i in range(max(0, n-29))]
    result = {"sequence_sha256": sequence_hash(seq), "length": n,
              "coordinate_system": "candidate_fragment_1_based_inclusive", "parameters": dict(PARAMETERS),
              "gravy": sum(KD[a] for a in seq)/n,
              "hydrophobic_segments": _segments([(i+1,v) for i,v in enumerate(kd) if v >= PARAMETERS['hydropathy_alert']], 19),
              "low_complexity_segments": _segments([(i+1,v) for i,v in enumerate(entropies) if v <= PARAMETERS['entropy_alert']], 12, minimum=True),
              "stp_enriched_segments": _segments([(i+1,v) for i,v in enumerate(stp) if v >= PARAMETERS['stp_alert']], 30),
              "algorithm": "Kyte-Doolittle sliding mean; Shannon window entropy (not SEG or disorder prediction)",
              "algorithm_version": "1", "functional_calibration": "not_performed"}
    result['biopython'] = {"status": "not_run"}
    if use_biopython:
        try:
            import Bio
            from Bio.SeqUtils.ProtParam import ProteinAnalysis
            from Bio.SeqUtils.ProtParamData import kd as bio_kd
            p = ProteinAnalysis(seq)
            profile = p.protein_scale(bio_kd, window=19) if n >= 19 else []
            agreement = abs(p.gravy()-result['gravy']) < 1e-8 and all(abs(a-b)<1e-8 for a,b in zip(kd,profile)) and len(kd)==len(profile)
            if not agreement:
                raise ValueError("Biopython hydropathy cross-check disagreed")
            result['biopython'] = {"status": "computed", "version": Bio.__version__, "kd_crosscheck": "pass",
                "molecular_weight": p.molecular_weight(), "isoelectric_point": p.isoelectric_point(),
                "charge_at_ph7": p.charge_at_pH(7), "source": "https://biopython.org/docs/latest/api/Bio.SeqUtils.ProtParam.html"}
        except ImportError:
            result['biopython'] = {"status": "not_available", "reason": "Optional biopython extra not installed"}
    return result


def validate_prediction(record, candidate, protein=None):
    """Accept normalized local outputs only, bound to this exact fragment.

The caller retains rejected records. No coordinates from a full protein are
silently transferred to a truncated fragment.
"""
    if record.get('software') not in {'IUPred2A', 'IUPred3', 'AIUPred', 'DeepTMHMM', 'DeepTMHMM2', 'SignalP'}:
        raise ValueError('Unsupported predictor')
    scope = record.get('input_scope')
    if record.get('coordinate_system') != '1-based-inclusive' or scope not in {'candidate_fragment', 'full_reference'}:
        raise ValueError('Predictor input scope/coordinates must explicitly match fragment or full reference')
    if scope == 'candidate_fragment':
        if record.get('candidate_id') != candidate.get('candidate_id') or record.get('sequence_sha256') != sequence_hash(candidate['sequence']):
            raise ValueError('Predictor identity/sequence hash mismatch')
        n = len(candidate['sequence'])
        start, end = 1, n
    else:
        from .common import interval, sourced
        if not protein or not sourced(protein):
            raise ValueError('Full-reference predictor requires sourced reference')
        ref = sequence(protein['sequence'])
        start, end = interval(candidate, len(ref))
        if record.get('accession') != protein.get('accession') or record.get('sequence_sha256') != sequence_hash(ref):
            raise ValueError('Full-reference predictor identity/sequence hash mismatch')
        if ref[start-1:end] != candidate['sequence'] or candidate.get('reference_sha256') != digest(ref):
            raise ValueError('Candidate/reference projection conflict')
        n = len(ref)
    if not record.get('version') or not isinstance(record.get('parameters'), dict) or not record.get('source'):
        raise ValueError('Predictor version, parameters and source required')
    checksum = record.get('raw_output_sha256', '')
    if len(checksum) != 64 or any(c not in '0123456789abcdef' for c in checksum):
        raise ValueError('Raw predictor output SHA256 required')
    if record['software'] in {'IUPred2A', 'IUPred3', 'AIUPred'}:
        scores = record.get('scores')
        if not isinstance(scores, list) or len(scores) != n or any(type(v) not in {int,float} or not math.isfinite(v) or not 0 <= v <= 1 for v in scores):
            raise ValueError('Disorder output needs one finite 0..1 score per residue')
    else:
        if not isinstance(record.get('regions'), list):
            raise ValueError('Topology predictor needs explicit regions list')
        for r in record['regions']:
            if r.get('kind') not in {'transmembrane','signal_peptide','cytoplasmic','inside','outside','extracellular','lumenal','periplasmic','reentrant','interfacial','transit_peptide'} or type(r.get('start')) is not int or type(r.get('end')) is not int or not 1 <= r['start'] <= r['end'] <= n:
                raise ValueError('Invalid predicted region')
    bound = dict(record, evidence_kind='prediction', functional_validation='not_performed', record_sha256=digest(record),
                 evaluated_candidate_id=candidate.get('candidate_id'), projection_interval=[start,end],
                 interpretation_limit='Native outside may be organelle lumen; isolated fragment SignalP cannot establish fusion secretion.')
    if scope == 'full_reference':
        if record.get('scores') is not None:
            bound['scores'] = record['scores'][start-1:end]
        if record.get('regions') is not None:
            bound['regions'] = [dict(r, start=max(start,r['start'])-start+1,
                end=min(end,r['end'])-start+1, reference_interval=[r['start'],r['end']])
                for r in record['regions'] if r['start'] <= end and start <= r['end']]
            bound['not_retained_regions'] = [r for r in record['regions'] if r['end'] < start or r['start'] > end]
    return bound
