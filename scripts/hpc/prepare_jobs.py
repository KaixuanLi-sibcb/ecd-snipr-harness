#!/usr/bin/env python3
"""Freeze exact reference jobs for scheduled tools; never substitute residues."""
import argparse
import hashlib
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ecd_snipr.common import file_hash, read_json, write_json


def prepare(reference_set, output, pilot_limit=None):
    source = Path(reference_set).resolve()
    output = Path(output).resolve()
    if output.exists() or output == source or source in output.parents:
        raise ValueError('Use a new manifest outside the immutable reference set')
    if pilot_limit is not None and pilot_limit <= 0:
        raise ValueError('Positive pilot limit required')
    definition = read_json(source/'set_definition.json')
    resolved, seen, unresolved = [], set(), []
    for entry in definition['entries']:
        if not entry.get('protein_file'):
            unresolved.append({k:entry.get(k) for k in ('input_id','accession','status','reason')})
            continue
        acc = entry.get('accession')
        if acc in seen:
            continue
        seen.add(acc)
        resolved.append(entry)
    selected = resolved[:pilot_limit] if pilot_limit else resolved
    jobs, exclusions = [], []
    for entry in selected:
        acc = entry['accession']
        try:
            if not isinstance(acc, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,80}', acc):
                raise ValueError('Unsafe accession identifier')
            path = (source/entry['protein_file']).resolve()
            if not path.is_relative_to(source) or file_hash(path) != entry['protein_sha256']:
                raise ValueError('Reference path/hash conflict')
            protein = read_json(path)
            seq = protein['sequence']
            if protein['accession'] != acc or not seq or seq != seq.strip():
                raise ValueError('Reference identity/sequence conflict')
            item = {'accession':acc, 'gene':entry.get('gene'), 'sequence':seq,
                    'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),
                    'protein_path':str(path), 'protein_sha256':file_hash(path)}
            invalid = sorted(set(seq)-set('ACDEFGHIKLMNPQRSTVWY'))
            if invalid:
                exclusions.append(dict(item, status='unsupported_sequence', invalid_residues=invalid,
                                       reason='Strict 20-AA workflow contract; no substitution'))
            else:
                jobs.append(item)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            exclusions.append({'accession':acc, 'gene':entry.get('gene'),
                               'status':'reference_error', 'reason':str(exc)})
    manifest = {'reference_set_sha256':file_hash(source/'set_definition.json'),
                'reference_set':str(source), 'database':definition.get('database'),
                'source_completeness':definition.get('completeness','unknown'),
                'source_reference_total':len(resolved), 'source_input_rows':len(definition['entries']),
                'unresolved_source_rows':unresolved, 'pilot':pilot_limit is not None,
                'scope_total':len(selected), 'eligible_total':len(jobs),
                'excluded_total':len(exclusions), 'jobs':jobs, 'exclusions':exclusions}
    write_json(output, manifest)
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-set', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--pilot-limit', type=int)
    args = p.parse_args()
    result = prepare(args.reference_set, args.output, args.pilot_limit)
    print({k:result[k] for k in ('scope_total','eligible_total','excluded_total','pilot')})


if __name__ == '__main__':
    main()
