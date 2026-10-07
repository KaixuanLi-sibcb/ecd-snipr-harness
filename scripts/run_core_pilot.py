#!/usr/bin/env python3
"""Small public topology pilot; never a functional accuracy benchmark."""
import argparse
import json
from pathlib import Path

from ecd_snipr.acquisition import build_set, read_target_list
from ecd_snipr.common import file_hash, read_json, write_json
from ecd_snipr.harness import verify_bundle
from ecd_snipr.interpro import fetch_domains
from ecd_snipr.receiver_risk import audit_risk
from ecd_snipr.screening import run_screening


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', required=True, help='New directory outside the repository; all caches/results stay local')
    parser.add_argument('--live', action='store_true', help='Retrieve public accession annotations only; no private sequences submitted')
    args = parser.parse_args()
    root = Path(args.outdir).resolve()
    repository = Path(__file__).resolve().parents[1]
    if root == repository or repository in root.parents:
        raise ValueError('Pilot must be outside the repository')
    root.mkdir(parents=True, exist_ok=False)
    definition = build_set(rows=read_target_list(repository/'examples/public_core_pilot.tsv'),
                           cache=root/'uniprot-cache', outdir=root/'analysis-set', offline=not args.live)
    screened = run_screening(root/'analysis-set',root/'screening',pilot=True)
    run = Path(screened['run_dir'])
    frozen = file_hash(run/'manifest.json')
    # Include no-candidate references in source coverage, not only easy successes.
    source_records = {}
    topologies = {}
    for entry in definition['entries']:
        if entry.get('protein_file'):
            p = read_json(root/'analysis-set'/entry['protein_file'])
            topologies[p['accession']] = p['topology']
            r = fetch_domains(p['accession'],root/'interpro-cache',offline=not args.live)
            source_records[p['accession']] = {k:v for k,v in r.items() if k != 'snapshot'}
    result = audit_risk(run,root/'analysis-set',root/'core-evidence',interpro_cache=root/'interpro-cache')
    overlay = Path(result['run_dir'])
    checks = {'screen_manifest':verify_bundle(run), 'overlay_manifest':verify_bundle(overlay),
              'source_unchanged':file_hash(run/'manifest.json') == frozen,
              'design_content_unchanged':[(c['candidate_id'],c.get('sequence'),c.get('screening_recommendation')) for c in read_json(run/'candidates.json')] ==
                                        [(c['candidate_id'],c.get('sequence'),c.get('screening_recommendation')) for c in read_json(overlay/'candidates.json')]}
    complete = definition['completeness'] == 'complete' and result['summary']['completeness'] == 'complete' and all(r['status'] in {'cached','fetched'} for r in source_records.values())
    summary = {'scope':'public_four_reference_topology_pilot_not_full_cohort', 'counts':screened['summary']['counts'],
               'topologies':topologies,'interpro_sources':source_records,'core_summary':result['summary']['core_evidence'],
               'checks':checks,'completeness':'complete' if complete else 'partial',
               'prediction_execution':{name:'not_run' for name in ('SignalP','DeepTMHMM','IUPred')},
               'functional_accuracy':'not_estimated','screen_run':str(run),'core_run':str(overlay)}
    write_json(root/'pilot_summary.json',summary)
    report = ['# Public core-evidence pilot', '', '**Software/annotation validation only, not biological accuracy.**', '',
              f"References resolved: {summary['counts']['resolved_identities']}/4. Completeness: {summary['completeness']}.",
              f"Usable candidates: {summary['counts']['candidates_usable']}. Design content unchanged: {checks['design_content_unchanged']}.", '',
              'Topology classes: '+json.dumps(topologies,sort_keys=True), '',
              'InterPro source coverage includes no-candidate inputs. Optional predictors were not executed.',
              'See pilot_summary.json and candidate_domain_coverage.tsv for exact source versions, sequence hashes and intervals.']
    (root/'PILOT_REPORT.md').write_text('\n'.join(report)+'\n')
    print(json.dumps(summary,ensure_ascii=False,indent=2))
    return 0 if complete and all(checks.values()) else 3


if __name__ == '__main__':
    raise SystemExit(main())
