#!/usr/bin/env python3
"""Scheduled local API inference with immutable records and exact-sequence import."""
import argparse
from collections import Counter
import gc
import importlib.metadata
import json
import os
from pathlib import Path
import signal
import socket
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ecd_snipr.common import file_hash, now, write_json
from ecd_snipr.local_predictors import _check_install
from ecd_snipr.predictor_import import import_prediction

def timeout_handler(signum, frame):
    raise TimeoutError('Per-reference inference deadline exceeded')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs', required=True, type=Path)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--outdir', required=True, type=Path)
    parser.add_argument('--tool', choices=['IUPred2A','DeepTMHMM2'], required=True)
    parser.add_argument('--shard', type=int, default=1)
    parser.add_argument('--shards', type=int, default=1)
    parser.add_argument('--device', choices=['cpu','cuda'], default='cpu')
    parser.add_argument('--per-reference-timeout', type=int, default=900)
    args = parser.parse_args()
    if not 1 <= args.shard <= args.shards or args.per_reference_timeout <= 0:
        parser.error('Invalid shard/deadline')
    config = json.loads(args.config.read_text())
    tool = config['tools'][args.tool]
    if Path(sys.executable).absolute() != Path(tool['python']).absolute():
        raise ValueError('Run this worker with the configured tool Python')
    _check_install(tool)
    manifest = json.loads(args.jobs.read_text())
    jobs = sorted(manifest['jobs'],key=lambda j:j['accession'])[args.shard-1::args.shards]
    root = args.outdir.resolve()/args.tool
    root.mkdir(parents=True,exist_ok=True)
    contract = {'jobs_sha256':file_hash(args.jobs),'config_sha256':file_hash(args.config),
                'worker_sha256':file_hash(__file__),'tool':args.tool,'shard':args.shard,
                'shards':args.shards,'device':args.device,'per_reference_timeout':args.per_reference_timeout}
    parameters = {'execution':'scheduled_local_only','device':args.device,
                  'input_scope':'full_reference','simplify_io':True,
                  'adapter':'author_cli_single_batch_unpacking',
                  'adapter_sha256':file_hash(Path(__file__).with_name('dtm2_api_compat.py')),
                  'precision':'author_cli_autocast'} if args.tool=='DeepTMHMM2' else {
                  'execution':'scheduled_local_only','device':'cpu','input_scope':'full_reference',
                  'method':'iupred2','mode':'long'}
    if args.tool == 'DeepTMHMM2':
        contract['adapter_sha256'] = parameters['adapter_sha256']
    if args.tool == 'DeepTMHMM2':
        os.environ['TORCH_HOME'] = str(Path(tool['model_dir'])/'torch')
        import torch
        torch.set_num_threads(int(os.environ.get('NSLOTS','4')))
        if args.device == 'cuda' and not torch.cuda.is_available():
            raise RuntimeError('CUDA requested but unavailable; no silent CPU fallback')
        from dtm2_api_compat import ReferencePredictor
        predictor = ReferencePredictor(model_dir=tool['model_dir'],device=args.device)
    else:
        from iupred.iupred2a import iupred
    states = []
    signal.signal(signal.SIGALRM,timeout_handler)
    for i, job in enumerate(jobs,1):
        acc = job['accession']
        folder = root/acc
        folder.mkdir(exist_ok=True)
        status_file = folder/'status.json'
        started = time.monotonic()
        try:
            protein = json.loads(Path(job['protein_path']).read_text())
            if file_hash(job['protein_path']) != job['protein_sha256'] or protein['sequence'] != job['sequence']:
                raise ValueError('Source protein changed')
            if status_file.exists():
                prior = json.loads(status_file.read_text())
                if prior.get('contract') != contract:
                    raise ValueError('Resume contract changed; use a new output directory')
                if prior['status'] == 'computed':
                    record = json.loads((folder/'evidence.json').read_text())
                    reproduced = import_prediction(args.tool,record['raw_output_path'],record['raw_input_path'],
                                                   protein,tool['version'],record['parameters'])
                    if any(record[k]!=v for k,v in reproduced.items()):
                        raise ValueError('Resume evidence/hash conflict')
                    states.append(prior)
                    continue
                # Failed attempts are preserved; retries require a new output root.
                states.append(prior)
                continue
            fasta = folder/'input.fasta'
            if fasta.exists():
                if fasta.read_text() != '>'+acc+'\n'+job['sequence']+'\n':
                    raise ValueError('Input FASTA conflict')
            else:
                fasta.write_text('>'+acc+'\n'+job['sequence']+'\n')
            signal.alarm(args.per_reference_timeout)
            if args.tool == 'DeepTMHMM2':
                predicted = predictor.predict({acc:job['sequence']},simplify_io=True,progress_bar=False)
                raw = folder/'predicted_topologies.3line'
                predicted.to_3line(raw)
                predicted.to_json(folder/'vendor_prediction.json')
            else:
                scores,_ = iupred(job['sequence'],mode='long')
                if len(scores) != len(job['sequence']):
                    raise ValueError('Incomplete IUPred score vector')
                raw = folder/'scores.tsv'
                raw.write_text('# IUPred2 API output: '+acc+'\n# pos aa iupred2\n'+''.join(
                      f'{p}\t{aa}\t{float(score):.17g}\n' for p,(aa,score) in enumerate(zip(job['sequence'],scores),1)))
            signal.alarm(0)
            record = import_prediction(args.tool,raw,fasta,protein,tool['version'],parameters)
            record['execution_provenance'] = {'configuration_sha256':file_hash(args.config),
                 'code_data_sha256':tool['code_data_sha256'],'model_sha256':tool['models'],
                 'worker_sha256':contract['worker_sha256'],'job_id':os.environ.get('JOB_ID'),
                 'host':socket.gethostname(),'executed_at':now(),
                 'python':tool['python'], 'executable_sha256':tool['executable_sha256']}
            write_json(folder/'evidence.json',record)
            state = {'accession':acc,'status':'computed','contract':contract,'seconds':time.monotonic()-started}
        except Exception as exc:
            signal.alarm(0)
            if status_file.exists():
                raise
            state = {'accession':acc,'status':'error','reason':type(exc).__name__+': '+str(exc),
                     'contract':contract,'seconds':time.monotonic()-started}
            if args.tool == 'DeepTMHMM2':
                gc.collect()
                if args.device == 'cuda':
                    torch.cuda.empty_cache()
        write_json(status_file,state)
        states.append(state)
        print(json.dumps(dict(state,index=i,total=len(jobs))),flush=True)
    _check_install(tool)
    if args.tool == 'DeepTMHMM2' and file_hash(Path(__file__).with_name('dtm2_api_compat.py')) != contract['adapter_sha256']:
        raise ValueError('Compatibility adapter changed during inference')
    summary = {'scope':'pilot' if manifest.get('pilot') else 'declared_reference_set',
               'software':args.tool,'states':dict(Counter(s['status'] for s in states)),
               'shard':args.shard,'shards':args.shards,'total':len(jobs),'contract':contract,
               'functional_accuracy':'not_estimated'}
    write_json(root/f'shard-{args.shard:02d}-summary.json',summary)
    print(json.dumps(summary),flush=True)
    return int(any(s['status']!='computed' for s in states))

if __name__ == '__main__':
    raise SystemExit(main())
