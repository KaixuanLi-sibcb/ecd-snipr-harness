"""Read local vendor outputs; no model installation, network or guessed results."""
from pathlib import Path
import re

from .common import file_hash, sequence
from .sequence_tools import sequence_hash


def import_prediction(software, raw_path, input_fasta, protein, version, parameters):
    accession = protein.get('accession')
    if not accession or not version or not isinstance(parameters, dict):
        raise ValueError('Accession, actual software version and parameters required')
    fasta = Path(input_fasta).read_text().splitlines()
    lines = [line.strip() for line in fasta if line.strip()]
    if not lines or not lines[0].startswith('>') or sum(line.startswith('>') for line in lines) != 1:
        raise ValueError('One explicit reference sequence per import is required')
    if lines[0][1:].split()[0] != accession:
        raise ValueError('FASTA identifier must be exact reference accession')
    ref = sequence(''.join(lines[1:]))
    if ref != sequence(protein['sequence']):
        raise ValueError('Predictor input FASTA/reference sequence conflict')
    text = Path(raw_path).read_text()
    record = {'software': software, 'version': version, 'parameters': parameters,
              'accession': accession, 'sequence_sha256': sequence_hash(ref),
              'input_scope': 'full_reference', 'coordinate_system': '1-based-inclusive',
              'raw_output_path': str(Path(raw_path).resolve()), 'raw_output_sha256': file_hash(raw_path),
              'raw_input_path': str(Path(input_fasta).resolve()), 'raw_input_sha256': file_hash(input_fasta)}
    if software in {'IUPred2A', 'IUPred3', 'AIUPred'}:
        record['source'] = 'https://iupred3.elte.hu/'
        residues, scores = [], []
        for line in text.splitlines():
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 3 or int(parts[0]) != len(scores) + 1:
                raise ValueError('IUPred position/order conflict')
            residues.append(parts[1])
            scores.append(float(parts[2]))
        if ''.join(residues) != ref:
            raise ValueError('IUPred residue/reference conflict')
        record['scores'] = scores
    elif software in {'DeepTMHMM', 'DeepTMHMM2'}:
        record['source'] = ('https://github.com/fteufel/DeepTMHMM2' if software == 'DeepTMHMM2'
                            else 'https://dtu.biolib.com/DeepTMHMM')
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        if len(lines) != 3 or not lines[0].startswith('>') or lines[0][1:].split()[0] != accession or lines[1] != ref or len(lines[2]) != len(ref):
            raise ValueError('Expected one DeepTMHMM three-line record with exact sequence')
        labels = lines[2]
        mapping = {'M':'transmembrane', 'B':'transmembrane', 'S':'signal_peptide', 'I':'inside', 'O':'outside', 'P':'periplasmic'}
        if software == 'DeepTMHMM2':
            if parameters.get('simplify_io') is not True:
                raise ValueError('DeepTMHMM2 import requires explicit simplify_io; membrane-specific labels are not guessed')
            mapping = dict(mapping, R='reentrant', F='interfacial')
            mapping['>'] = 'transit_peptide'
            del mapping['P']
        if set(labels) - set(mapping):
            raise ValueError('Unsupported topology labels; no silent DeepTMHMM2 format assumption')
        regions = []
        for i, label in enumerate(labels, 1):
            if regions and regions[-1]['raw_label'] == label:
                regions[-1]['end'] = i
            else:
                regions.append({'kind': mapping[label], 'start': i, 'end': i, 'raw_label': label})
        record['regions'] = regions
    elif software == 'SignalP':
        record['source'] = 'https://services.healthtech.dtu.dk/services/SignalP-6.0/'
        if any('Prediction' in line and 'CS Position' in line for line in text.splitlines() if line.startswith('#')):
            lines = [line for line in text.splitlines() if line.strip() and not line.startswith('#')]
            if len(lines) != 1:
                raise ValueError('One SignalP summary record per exact FASTA required')
            parts = lines[0].split('\t')
            if len(parts) < 3 or parts[0].strip() != accession:
                raise ValueError('SignalP summary identity/format conflict')
            prediction = parts[1].strip()
            record['predicted_class'] = prediction
            if prediction == 'OTHER':
                if re.search(r'CS pos', lines[0], re.I):
                    raise ValueError('No-SP class contradicts reported cleavage')
                record['regions'] = []
                record['signal_peptide_state'] = 'predicted_absent_not_experimentally_absent'
            elif prediction == 'SP':
                match = re.search(r'CS pos:\s*(\d+)\s*-\s*(\d+)', lines[0])
                if not match or int(match[2]) != int(match[1])+1:
                    raise ValueError('SP class requires an explicit adjacent cleavage site')
                record['regions'] = [{'kind':'signal_peptide','start':1,'end':int(match[1])}]
                record['signal_peptide_state'] = 'predicted_present'
            else:
                raise ValueError('Unsupported SignalP class; human eukarya SP/OTHER only')
            from .sequence_tools import validate_prediction
            from .common import digest
            validate_prediction(record,{'candidate_id':'import_validation_only','sequence':ref,'start':1,
                'end':len(ref),'reference_sha256':digest(ref)},protein)
            return record
        regions = []
        saw_id = False
        for line in text.splitlines():
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) != 9 or parts[0] != accession:
                raise ValueError('SignalP GFF3 identity/format conflict')
            saw_id = True
            if parts[2].lower() in {'signal_peptide', 'signalpeptide'}:
                regions.append({'kind':'signal_peptide', 'start':int(parts[3]), 'end':int(parts[4])})
        # Empty GFF is not a checked negative: import positive regions only.
        if not saw_id or not regions:
            raise ValueError('No explicit signal region in GFF; a no-SP summary needs a separately supported format')
        record['regions'] = regions
    else:
        raise ValueError('Unsupported predictor')
    from .sequence_tools import validate_prediction
    from .common import digest
    dummy = {'candidate_id':'import_validation_only', 'sequence':ref, 'start':1, 'end':len(ref), 'reference_sha256':digest(ref)}
    validate_prediction(record,dummy,protein)
    return record
