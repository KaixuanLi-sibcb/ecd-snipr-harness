#!/usr/bin/env python3
"""Independently reimport every raw prediction; coverage is not biological accuracy."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ecd_snipr.common import file_hash, now, tsv, write_json
from ecd_snipr.predictor_import import import_prediction

def verify_record(record, protein):
    provenance=record.get('execution_provenance',{})
    for part in ('batch_input','batch_output'):
        if part+'_path' in provenance and file_hash(provenance[part+'_path'])!=provenance[part+'_sha256']:
            raise ValueError('Batch source changed')
    if 'batch_output_path' in provenance:
        lines=Path(provenance['batch_output_path']).read_text().splitlines()
        header='\n'.join(line for line in lines if line.startswith('#'))+'\n'
        selected=[line for line in lines if line.strip() and not line.startswith('#')
                  and line.split('\t')[0].strip()==record['accession']]
        if len(selected)!=1 or Path(record['raw_output_path']).read_text()!=header+selected[0]+'\n':
            raise ValueError('Separated row differs from original vendor batch')
    reproduced = import_prediction(record['software'],record['raw_output_path'],record['raw_input_path'],
                                   protein,record['version'],record['parameters'])
    if any(record.get(k)!=v for k,v in reproduced.items()):
        raise ValueError('Independent reimport differs')
    return reproduced

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs',required=True,type=Path)
    parser.add_argument('--run',required=True,type=Path)
    parser.add_argument('--outdir',required=True,type=Path)
    parser.add_argument('--tools',nargs='+',default=['IUPred2A','DeepTMHMM2','SignalP'])
    parser.add_argument('--tool-root',action='append',default=[],help='Explicit TOOL=/path/to/tool-results override')
    args = parser.parse_args()
    if not args.tools or len(set(args.tools))!=len(args.tools) or set(args.tools)-{'IUPred2A','DeepTMHMM2','SignalP'}:
        parser.error('Unique supported tools required')
    args.outdir.mkdir(parents=True,exist_ok=False)
    manifest = json.loads(args.jobs.read_text())
    source = Path(manifest['reference_set'])
    if file_hash(source/'set_definition.json') != manifest['reference_set_sha256']:
        raise ValueError('Frozen reference set changed')
    entries = manifest['jobs']+manifest['exclusions']
    if len({j['accession'] for j in entries})!=len(entries):
        raise ValueError('Duplicate reference in jobs manifest')
    rows,records = [],[]
    overrides = dict(item.split('=',1) for item in args.tool_root)
    if set(overrides)-set(args.tools):
        parser.error('Unknown tool root override')
    for tool in args.tools:
        for job in entries:
            acc=job['accession']
            row={'accession':acc,'gene':job.get('gene'),'software':tool,'status':'not_run','reason':''}
            if job.get('status'):
                rows.append(dict(row,status=job['status'],reason=job['reason']))
                continue
            folder=Path(overrides.get(tool,args.run/tool))/acc
            status=folder/'status.json'
            if not status.exists():
                rows.append(dict(row,reason='No per-reference completion record'))
                continue
            state=json.loads(status.read_text())
            if state['accession']!=acc:
                raise ValueError('Status accession conflict')
            if state['status']!='computed':
                rows.append(dict(row,status=state['status'],reason=state.get('reason','')))
                continue
            try:
                if state.get('contract',{}).get('jobs_sha256')!=file_hash(args.jobs):
                    raise ValueError('Completion record belongs to different frozen jobs')
                if file_hash(job['protein_path'])!=job['protein_sha256']:
                    raise ValueError('Protein source changed')
                protein=json.loads(Path(job['protein_path']).read_text())
                record=json.loads((folder/'evidence.json').read_text())
                if record['software']!=tool or record['accession']!=acc:
                    raise ValueError('Record identity conflict')
                verify_record(record,protein)
                records.append(record)
                rows.append(dict(row,status='computed',raw_output_sha256=record['raw_output_sha256']))
            except (OSError,ValueError,KeyError,TypeError) as exc:
                rows.append(dict(row,status='integrity_error',reason=str(exc)))
    states={tool:dict(Counter(r['status'] for r in rows if r['software']==tool)) for tool in args.tools}
    eligible_complete=bool(manifest['eligible_total']) and all(states[t].get('computed',0)==manifest['eligible_total'] for t in args.tools)
    complete=eligible_complete and not manifest['excluded_total'] and not manifest.get('pilot') and not manifest.get('unresolved_source_rows') and manifest.get('source_completeness')=='complete'
    summary={'created_at':now(),'scope_denominator':manifest['scope_total'],
             'eligible_denominator':manifest['eligible_total'],
             'excluded_denominator':manifest['excluded_total'],
             'unsupported_denominator':sum(j.get('status')=='unsupported_sequence' for j in manifest['exclusions']),
             'source_completeness':manifest.get('source_completeness','unknown'),
             'source_reference_denominator':manifest.get('source_reference_total',manifest['scope_total']),
             'pilot':bool(manifest.get('pilot')),
             'requested_tools':args.tools,'states':states,
             'execution_complete_for_eligible':eligible_complete,
             'completeness':'complete' if complete else 'partial',
             'functional_accuracy':'not_estimated','experimental_labels_used':False,
             'validation':'Exact reference, source hash, raw output and parser reproduction; not functional SNIPR validation',
             'jobs_sha256':file_hash(args.jobs),'verifier_sha256':file_hash(__file__)}
    write_json(args.outdir/'tool_evidence.json',records)
    write_json(args.outdir/'summary.json',summary)
    write_json(args.outdir/'predictor_summary.json',summary)
    tsv(args.outdir/'predictor_status.tsv',rows,['accession','gene','software','status','reason','raw_output_sha256'])
    write_json(args.outdir/'manifest.json',{'state':'complete','artifacts':{
              p.name:file_hash(p) for p in args.outdir.iterdir() if p.is_file()},
              'nested_artifacts':{},'analysis_completeness':summary['completeness'],
              'inputs':dict(summary,reference_set_path=str(source),
                            set_sha256=manifest['reference_set_sha256'])})
    print(json.dumps(summary,indent=2),flush=True)
    return int(not eligible_complete or any(s.get('integrity_error',0) for s in states.values()))

if __name__=='__main__':
    raise SystemExit(main())
