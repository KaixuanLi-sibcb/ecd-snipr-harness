"""Exact-fragment boundary/domain cross-checks; never a receptor predictor."""
from collections import Counter
from pathlib import Path

from .common import digest, interval, overlaps, sourced, sequence, tsv, write_json
from .interpro import normalized_domains
from .sequence_tools import sequence_hash, validate_prediction

CORE_POLICY = {'version': '1', 'coordinates': '1-based-inclusive',
               'family_is_domain': False, 'prediction_overwrites_boundary': False,
               'missing_is_negative': False, 'functional_accuracy': 'not_established'}


def check_core(candidate, protein=None, interpro=None, predictions=()):
    result = {'policy': CORE_POLICY, 'candidate_id': candidate.get('candidate_id'),
              'status': 'not_assessed', 'checks': [], 'domains': [], 'topology_predictions': [],
              'rejected_evidence': [], 'reason_codes': [],
              'coverage': {'uniprot': 'not_supplied', 'interpro': 'not_run',
                           'SignalP': 'not_run', 'DeepTMHMM': 'not_run', 'disorder': 'not_run'},
              'claim_limit': 'Boundary and structural-unit checks do not establish surface expression, epitope preservation or SNIPR function.'}
    if not protein or candidate.get('design_status') == 'blocked':
        return result
    try:
        reference, seq = sequence(protein['sequence']), sequence(candidate['sequence'])
        start, end = interval(candidate, len(reference))
        if reference[start-1:end] != seq or candidate.get('reference_sha256') != digest(reference) or not sourced(protein):
            raise ValueError('reference_sequence_or_provenance_conflict')
    except (ValueError, KeyError, TypeError) as exc:
        return dict(result, status='reference_conflict', reason_codes=[str(exc)])
    result['reference_sequence_sha256'] = sequence_hash(reference)
    result['candidate_sequence_sha256'] = sequence_hash(seq)
    result['coverage']['uniprot'] = 'evaluated'
    result['checks'].append({'code': 'exact_reference_fragment', 'state': 'confirmed_sequence_correspondence'})
    domains, extras = [], []
    for index, feature in enumerate(protein.get('features', [])):
        if feature.get('kind') not in {'domain', 'repeat', 'transmembrane', 'intramembrane', 'signal_peptide', 'cytoplasmic'}:
            continue
        try:
            a, b = interval(feature, len(reference))
            if not sourced(feature):
                raise ValueError('unsourced_feature')
        except ValueError as exc:
            result['rejected_evidence'].append({'feature_index': index, 'reason': str(exc)})
            continue
        if feature['kind'] in {'domain', 'repeat'}:
            domains.append(dict(feature, database='UniProt', fragments=[{'start': a, 'end': b}]))
        elif overlaps((start, end), (a, b)):
            result['checks'].append({'code': 'retained_native_exclusion_region', 'state': 'conflict',
                                    'kind': feature['kind'], 'reference_interval': [a, b], 'evidence': feature['evidence']})
    if interpro:
        result['coverage']['interpro'] = interpro['status']
        if interpro['status'] in {'cached', 'fetched'}:
            try:
                extras = normalized_domains(interpro, protein)
                domains.extend(extras)
                result['coverage']['interpro'] = 'evaluated'
                result['checks'].append({'code': 'interpro_sequence_correspondence', 'state': 'confirmed_sequence_correspondence',
                                        'raw_response_sha256': interpro['raw_response_sha256'], 'version': interpro['version']})
            except (ValueError, KeyError, TypeError) as exc:
                result['coverage']['interpro'] = 'conflict'
                result['rejected_evidence'].append({'source': 'InterPro', 'reason': str(exc)})
        else:
            result['rejected_evidence'].append({'source': 'InterPro', 'reason': interpro.get('reason', 'source_unavailable')})
    for domain in domains:
        fragments = domain['fragments']
        if all(start <= f['start'] <= f['end'] <= end for f in fragments):
            state = 'retained'
        elif any(overlaps((start, end), (f['start'], f['end'])) for f in fragments):
            state = 'cut'
        else:
            state = 'omitted'
        result['domains'].append(dict(domain, retention=state,
            autonomous_domain_not_established=domain.get('kind') in {'repeat', 'homologous_superfamily'},
            discontinuous=len(fragments) > 1 or any(f.get('dc-status', 'CONTINUOUS') != 'CONTINUOUS' for f in fragments)))
        if state == 'cut':
            result['checks'].append({'code': 'domain_boundary_cut', 'state': 'review',
                                    'database': domain['database'], 'name': domain.get('name'),
                                    'reference_interval': [domain['start'], domain['end']], 'evidence': domain['evidence']})
    # Repeated InterPro/member hits stay inspectable; they are not independent votes.
    result['checks'].append({'code': 'domain_annotations_available', 'state': 'evaluated' if domains else 'missing',
                            'annotation_records': len(domains), 'independent_vote_count': 'not_used'})
    for prediction in predictions:
        try:
            bound = validate_prediction(prediction, candidate, protein)
            source = 'DeepTMHMM' if bound['software']=='DeepTMHMM2' else bound['software']
            result['coverage'][source if source in {'SignalP', 'DeepTMHMM'} else 'disorder'] = 'evaluated'
            if source in {'SignalP', 'DeepTMHMM'}:
                result['topology_predictions'].append(bound)
                for region in bound['regions']:
                    if region['kind'] in {'transmembrane', 'signal_peptide', 'cytoplasmic', 'inside','reentrant','interfacial','transit_peptide'}:
                        result['checks'].append({'code': 'predicted_exclusion_region_retained', 'state': 'review',
                                                'region': region, 'evidence': bound})
        except (ValueError, KeyError, TypeError) as exc:
            result['rejected_evidence'].append({'source': prediction.get('software'), 'reason': str(exc)})
            source = prediction.get('software')
            if source == 'DeepTMHMM2':
                source = 'DeepTMHMM'
            if source in {'SignalP','DeepTMHMM'}:
                result['coverage'][source] = 'rejected'
            elif source in {'IUPred2A','IUPred3','AIUPred'}:
                result['coverage']['disorder'] = 'rejected'
    result['reason_codes'] = sorted({c['code'] for c in result['checks'] if c['state'] in {'review', 'conflict'}})
    result['status'] = 'crosscheck_review_required' if result['reason_codes'] else 'checks_completed_no_conflict'
    if result['coverage']['interpro'] == 'conflict':
        result['status'] = 'crosscheck_source_conflict'
    if result['rejected_evidence'] and result['status'] == 'checks_completed_no_conflict':
        result['status'] = 'checks_partial_evidence_unresolved'
    return result


def export_core(outdir, candidates):
    root = Path(outdir)
    rows, domains = [], []
    for c in candidates:
        r = c.get('receiver_risk', {}).get('core_evidence') or check_core(c)
        rows.append({'candidate_id': c['candidate_id'], 'accession': c.get('accession'),
                     'start': c.get('start'), 'end': c.get('end'), 'status': r['status'],
                     'reason_codes': r['reason_codes'], 'coverage': r['coverage'],
                     'checks': r['checks'], 'rejected_evidence': r['rejected_evidence']})
        domains.extend(dict(candidate_id=c['candidate_id'], **d) for d in r['domains'])
    tsv(root / 'candidate_core_evidence.tsv', rows, ['candidate_id', 'accession', 'start', 'end', 'status', 'reason_codes', 'coverage', 'checks', 'rejected_evidence'])
    tsv(root / 'candidate_domain_coverage.tsv', domains, ['candidate_id', 'database', 'entry_accession', 'kind', 'name', 'start', 'end', 'fragments', 'retention', 'discontinuous', 'evidence'])
    summary = {'candidate_records': len(rows), 'statuses': dict(Counter(r['status'] for r in rows)),
               'coverage_states': {source: dict(Counter(r['coverage'][source] for r in rows)) for source in ('uniprot', 'interpro', 'SignalP', 'DeepTMHMM', 'disorder')},
               'functional_accuracy': 'not_established', 'policy': CORE_POLICY}
    write_json(root / 'core_evidence_summary.json', summary)
    return summary
