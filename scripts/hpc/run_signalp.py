#!/usr/bin/env python3
"""Run licensed, unmodified SignalP6 in CPU batches with auditable row separation."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import socket
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ecd_snipr.common import file_hash,now,write_json
from ecd_snipr.local_predictors import _check_install
from ecd_snipr.predictor_import import import_prediction

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs',required=True,type=Path)
    parser.add_argument('--config',required=True,type=Path)
    parser.add_argument('--outdir',required=True,type=Path)
    parser.add_argument('--batch-size',type=int,default=512)
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--batch-timeout',type=int,default=7200)
    args=parser.parse_args()
    if min(args.batch_size,args.threads,args.batch_timeout)<=0:
        parser.error('Positive batch size, thread count and timeout required')
    tool=json.loads(args.config.read_text())['tools']['SignalP']
    _check_install(tool)
    manifest=json.loads(args.jobs.read_text())
    jobs=sorted(manifest['jobs'],key=lambda j:j['accession'])
    root=args.outdir.resolve()/'SignalP'
    root.mkdir(parents=True,exist_ok=True)
    parameters={'execution':'scheduled_local_only','device':'cpu','input_scope':'full_reference',
                'mode':'fast','organism':'eukarya','format':'none','write_procs':1,'torch_num_threads':args.threads}
    contract={'jobs_sha256':file_hash(args.jobs),'config_sha256':file_hash(args.config),
              'worker_sha256':file_hash(__file__),'batch_size':args.batch_size,
              'batch_timeout':args.batch_timeout,'parameters':parameters}
    states=[]
    for start in range(0,len(jobs),args.batch_size):
        selected=jobs[start:start+args.batch_size]
        batch=root/f'batch-{start//args.batch_size+1:03d}'
        batch.mkdir(exist_ok=True)
        completion=batch/'complete.json'
        if (batch/'error.json').exists():
            raise ValueError('Preserve failed batch; select a new output root')
        if completion.exists():
            previous=json.loads(completion.read_text())
            if previous['contract']!=contract:
                raise ValueError('Resume contract conflict; select a new output root')
            if file_hash(batch/'input.fasta')!=previous['input_sha256'] or file_hash(batch/'vendor/prediction_results.txt')!=previous['output_sha256']:
                raise ValueError('Resume batch hash conflict; preserve existing results')
        try:
            if completion.exists():
                previous=json.loads(completion.read_text())
                if previous['contract']!=contract: raise ValueError('Resume contract conflict')
                if file_hash(batch/'input.fasta')!=previous['input_sha256'] or file_hash(batch/'vendor/prediction_results.txt')!=previous['output_sha256']:
                    raise ValueError('Resume batch hash conflict')
            elif (batch/'invocation.json').exists():
                raise ValueError('Incomplete prior batch; preserve and retry under new output root')
            else:
                fasta=batch/'input.fasta'
                fasta.write_text(''.join('>'+j['accession']+'\n'+j['sequence']+'\n' for j in selected))
                command=[tool['executable'],'--fastafile',str(fasta),'--output_dir',str(batch/'vendor'),
                         '--model_dir',tool['model_dir'],'--organism','eukarya','--mode','fast',
                         '--format','none','--write_procs','1','--torch_num_threads',str(args.threads)]
                write_json(batch/'invocation.json',{'command':command,'contract':contract,'started_at':now()})
                with (batch/'stdout.log').open('w') as out,(batch/'stderr.log').open('w') as err:
                    proc=subprocess.run(command,stdout=out,stderr=err,timeout=args.batch_timeout)
                if proc.returncode: raise ValueError('SignalP exit '+str(proc.returncode))
                write_json(completion,{'contract':contract,'input_sha256':file_hash(fasta),
                          'output_sha256':file_hash(batch/'vendor/prediction_results.txt'),'completed_at':now()})
            raw_batch=batch/'vendor/prediction_results.txt'
            lines=raw_batch.read_text().splitlines()
            header='\n'.join(line for line in lines if line.startswith('#'))+'\n'
            data=[line for line in lines if line.strip() and not line.startswith('#')]
            lookup={}
            for line in data:
                acc=line.split('\t')[0].strip()
                if acc in lookup: raise ValueError('Duplicate vendor output identity')
                lookup[acc]=line
            if set(lookup)!={j['accession'] for j in selected}:
                raise ValueError('Batch identity coverage mismatch')
            batch_provenance={
                'batch_input_path':str(batch/'input.fasta'),'batch_input_sha256':file_hash(batch/'input.fasta'),
                'batch_output_path':str(raw_batch),'batch_output_sha256':file_hash(raw_batch)}
        except Exception as exc:
            write_json(batch/'error.json',{'reason':str(exc),'contract':contract,'timestamp':now()})
            for job in selected:
                folder=root/job['accession'];folder.mkdir(exist_ok=True)
                state={'accession':job['accession'],'status':'error','reason':str(exc),'contract':contract}
                write_json(folder/'status.json',state);states.append(state)
            print(json.dumps({'batch':str(batch),'status':'error','reason':str(exc)}),flush=True)
            continue
        for job in selected:
            acc=job['accession'];folder=root/acc;folder.mkdir(exist_ok=True)
            try:
                if file_hash(job['protein_path'])!=job['protein_sha256']: raise ValueError('Reference changed')
                protein=json.loads(Path(job['protein_path']).read_text())
                if protein['sequence']!=job['sequence'] or protein['accession']!=acc:
                    raise ValueError('Reference sequence/identity conflict')
                fasta=folder/'input.fasta';raw=folder/'prediction_results.txt'
                if not fasta.exists(): fasta.write_text('>'+acc+'\n'+job['sequence']+'\n')
                if not raw.exists(): raw.write_text(header+lookup[acc]+'\n')
                if fasta.read_text()!='>'+acc+'\n'+job['sequence']+'\n' or raw.read_text()!=header+lookup[acc]+'\n':
                    raise ValueError('Row separation conflict')
                record=import_prediction('SignalP',raw,fasta,protein,tool['version'],parameters)
                status_path=folder/'status.json'
                if status_path.exists():
                    prior=json.loads(status_path.read_text())
                    if prior.get('contract')!=contract or prior.get('status')!='computed':
                        raise ValueError('Preserve prior row; select a new output root')
                    existing=json.loads((folder/'evidence.json').read_text())
                    if any(existing.get(k)!=v for k,v in record.items()):
                        raise ValueError('Resume row evidence conflict')
                    states.append(prior)
                    continue
                record['execution_provenance']=dict(batch_provenance,configuration_sha256=file_hash(args.config),
                    code_data_sha256=tool['code_data_sha256'],model_sha256=tool['models'],
                    worker_sha256=contract['worker_sha256'],host=socket.gethostname(),job_id=os.environ.get('JOB_ID'),
                    executed_at=now(),executable_sha256=tool['executable_sha256'])
                write_json(folder/'evidence.json',record)
                state={'accession':acc,'status':'computed','contract':contract}
            except Exception as exc:
                if (folder/'status.json').exists():
                    raise
                state={'accession':acc,'status':'error','reason':str(exc),'contract':contract}
            write_json(folder/'status.json',state);states.append(state)
        print(json.dumps({'batch':str(batch),'processed':start+len(selected),'total':len(jobs),
              'states':dict(Counter(s['status'] for s in states))}),flush=True)
    _check_install(tool)
    summary={'software':'SignalP','total':len(jobs),'states':dict(Counter(s['status'] for s in states)),
             'scope':'pilot' if manifest.get('pilot') else 'declared_reference_set',
             'functional_accuracy':'not_estimated','contract':contract}
    write_json(root/'summary.json',summary)
    print(json.dumps(summary),flush=True)
    return int(any(s['status']!='computed' for s in states))

if __name__=='__main__':raise SystemExit(main())
