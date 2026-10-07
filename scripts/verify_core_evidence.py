#!/usr/bin/env python3
"""Independent result cross-check: re-derive intervals without the risk engine."""
import argparse
import json
from pathlib import Path

from ecd_snipr.common import digest, file_hash, read_json
from ecd_snipr.harness import verify_bundle


def verify_core(run_dir, reference_set):
    root, refs = Path(run_dir), Path(reference_set)
    issues = []
    if not verify_bundle(root):
        return {'passed':False,'issues':['artifact_integrity_failure']}
    definition = read_json(refs/'set_definition.json')
    entries = {e.get('accession'):e for e in definition['entries'] if e.get('protein_file')}
    checked, domain_records = 0, 0
    for c in read_json(root/'candidates.json'):
        if c.get('design_status') == 'blocked' or not c.get('sequence'):
            continue
        cid = c['candidate_id']
        entry = entries.get(c.get('accession') or c.get('protein_id'))
        if not entry:
            issues.append(cid+':missing_reference'); continue
        path = (refs/entry['protein_file']).resolve()
        if not path.is_relative_to(refs.resolve()) or file_hash(path) != entry['protein_sha256']:
            issues.append(cid+':reference_hash_conflict'); continue
        p = read_json(path)
        s, e = c['start'], c['end']
        if type(s) is not int or type(e) is not int or not 1 <= s <= e <= len(p['sequence']) or p['sequence'][s-1:e] != c['sequence'] or c.get('reference_sha256') != digest(p['sequence']):
            issues.append(cid+':exact_fragment_conflict'); continue
        risk = c.get('receiver_risk', {})
        if risk.get('functional_confidence') != 'not_established' or risk.get('success_probability','absent') is not None:
            issues.append(cid+':unsupported_functional_claim')
        core = risk.get('core_evidence')
        if not core:
            issues.append(cid+':missing_core_evidence'); continue
        for d in core['domains']:
            fragments = d['fragments']
            if any(type(f.get('start')) is not int or type(f.get('end')) is not int or not 1 <= f['start'] <= f['end'] <= len(p['sequence']) for f in fragments):
                issues.append(cid+':domain_coordinate_conflict'); continue
            expected = 'retained' if all(s <= f['start'] and f['end'] <= e for f in fragments) else 'cut' if any(max(s,f['start']) <= min(e,f['end']) for f in fragments) else 'omitted'
            if d['retention'] != expected:
                issues.append(cid+':domain_retention_conflict')
            domain_records += 1
        checked += 1
    source_file = root/'core_source_records.json'
    if source_file.exists():
        for accession, record in read_json(source_file).items():
            if record.get('status') not in {'cached','fetched'}:
                continue
            path = Path(record['raw_snapshot_path'])
            if file_hash(path) != record['raw_response_sha256']:
                issues.append(accession+':snapshot_hash_conflict'); continue
            snapshot = read_json(path)
            for retrieval in snapshot['retrievals']:
                raw = path.parent/retrieval['raw_response_file']
                if raw.parent != path.parent or file_hash(raw) != retrieval['response_sha256']:
                    issues.append(accession+':raw_response_hash_conflict')
    return {'passed':not issues,'usable_candidates_checked':checked,'domain_records_checked':domain_records,
            'issues':issues,'scope':'independent_sequence_interval_retention_and_artifact_checks_not_biological_validation'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',required=True)
    parser.add_argument('--reference-set',required=True)
    args=parser.parse_args()
    result=verify_core(args.run_dir,args.reference_set)
    print(json.dumps(result,indent=2))
    return 0 if result['passed'] else 1


if __name__=='__main__':
    try:
        raise SystemExit(main())
    except (OSError,ValueError,KeyError,TypeError) as exc:
        print(json.dumps({'passed':False,'reason':str(exc)}))
        raise SystemExit(1)
