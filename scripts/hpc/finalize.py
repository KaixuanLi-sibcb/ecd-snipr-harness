#!/usr/bin/env python3
"""Verify scheduled raw results and optionally overlay frozen candidate evidence."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ecd_snipr.common import read_json, write_json
from ecd_snipr.harness import verify_bundle
from ecd_snipr.local_predictors import verify_predictions
from ecd_snipr.receiver_risk import audit_risk


def finalize(prediction_dir, screening_run=None, overlay_outdir=None):
    root = Path(prediction_dir).resolve()
    if not verify_predictions(root):
        raise ValueError('Prediction bundle failed exact raw/reference verification')
    if bool(screening_run) != bool(overlay_outdir):
        raise ValueError('Supply both screening run and new overlay directory, or neither')
    result = {'prediction':read_json(root/'predictor_summary.json'),
              'functional_accuracy':'not_estimated'}
    result['completeness']=result['prediction']['completeness']
    if screening_run:
        if not verify_bundle(screening_run):
            raise ValueError('Frozen screening source failed verification')
        out = Path(overlay_outdir).resolve()
        source = Path(screening_run).resolve()
        setdir = Path(read_json(root/'manifest.json')['inputs']['reference_set_path'])
        if out.exists() or any(out==p or p in out.parents for p in (root,source,setdir)):
            raise ValueError('Select a new overlay root outside immutable inputs')
        candidates = read_json(source/'candidates.json')
        accessions = {c.get('accession') or c.get('protein_id') for c in candidates}
        records = read_json(root/'tool_evidence.json')
        selected = [r for r in records if r['accession'] in accessions]
        out.mkdir(parents=True)
        evidence = out/'candidate-tool-evidence.json'
        write_json(evidence, selected)
        overlay = audit_risk(source,setdir,out/'overlay',evidence)
        if not verify_bundle(overlay['run_dir']):
            raise ValueError('Candidate overlay failed artifact verification')
        result.update(overlay=overlay, candidate_reference_denominator=len(accessions),
                      records_for_candidate_references=len(selected),
                      records_without_candidate_reference=len(records)-len(selected))
        if overlay['summary']['completeness']!='complete':
            result['completeness']='partial'
        write_json(out/'acceptance.json',result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prediction-dir', type=Path, required=True)
    p.add_argument('--screening-run', type=Path)
    p.add_argument('--overlay-outdir', type=Path)
    args = p.parse_args()
    result = finalize(args.prediction_dir,args.screening_run,args.overlay_outdir)
    print(result)
    return 0 if result['completeness']=='complete' else 3


if __name__=='__main__':
    raise SystemExit(main())
