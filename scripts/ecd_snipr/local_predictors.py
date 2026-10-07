"""Local-only execution contract. No hosted inference or silent tool substitution."""
import json
import os
import re
import subprocess
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .common import digest, file_hash, now, read_json, sequence, tsv, write_json
from .predictor_import import import_prediction

TOOLS = ('IUPred2A', 'DeepTMHMM2', 'SignalP')


def verify_predictions(path):
    from .harness import verify_bundle
    from .receiver_risk import read_reference
    root = Path(path).resolve()
    if not verify_bundle(root):
        return False
    try:
        manifest = read_json(root/'manifest.json')
        source = Path(manifest['inputs']['reference_set_path'])
        if file_hash(source/'set_definition.json') != manifest['inputs']['set_sha256']:
            return False
        entries = {e['accession']:e for e in read_json(source/'set_definition.json')['entries'] if e.get('protein_file')}
        records = read_json(root/'tool_evidence.json')
        def load(acc):
            entry=entries[acc]; path=(source/entry['protein_file']).resolve()
            if not path.is_relative_to(source.resolve()):
                raise ValueError('Reference outside supplied set')
            return acc,read_reference(path,entry['protein_sha256'])
        with ThreadPoolExecutor(max_workers=4) as pool:
            proteins = dict(pool.map(load,dict.fromkeys(r['accession'] for r in records)))
        for name,checksum in manifest['nested_artifacts'].items():
            file = (root/name).resolve()
            if not file.is_relative_to(root) or file_hash(file)!=checksum:
                return False
        for record in records:
            for part in ('input','output'):
                if file_hash(record['raw_'+part+'_path']) != record['raw_'+part+'_sha256']:
                    return False
            reproduced = import_prediction(record['software'],record['raw_output_path'],record['raw_input_path'],
                proteins[record['accession']],record['version'],record['parameters'])
            if any(record.get(k)!=v for k,v in reproduced.items()):
                return False
        return True
    except (ValueError,KeyError,OSError,TypeError):
        return False


def _run(command, timeout=60, **kwargs):
    return subprocess.run(command, capture_output=True, text=True, timeout=timeout, **kwargs)


def inventory(python, package, model_dir=None):
    """Probe the installed package and hash code/data; does not load model weights."""
    command = [str(python), '-c',
        'import importlib.util,importlib.metadata,json;'
        f's=importlib.util.find_spec({package!r});'
        f'd={"signalp6" if package == "signalp" else package!r};'
        'print(json.dumps({"version":importlib.metadata.version(d),"distribution":d,"path":s.origin}))']
    result = _run(command)
    if result.returncode:
        raise ValueError('Package probe failed: '+result.stderr[-1000:])
    data = json.loads(result.stdout)
    root = Path(data['path']).resolve().parent
    files = {str(p.relative_to(root)): file_hash(p) for p in sorted(root.rglob('*'))
             if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pyc'}
    models = {}
    if model_dir:
        model = Path(model_dir).resolve()
        if not model.is_dir():
            raise ValueError('Model directory absent')
        models = {str(p.resolve()):file_hash(p) for p in sorted(model.rglob('*')) if p.is_file()}
        if not models:
            raise ValueError('Empty model directory')
    return dict(data, package=package, python=str(Path(python).absolute()),
                code_data_sha256=digest(files), files=files, models=models,
                model_dir=str(Path(model_dir).resolve()) if model_dir else None)


def configure(output, iupred_python=None, dtm_python=None, dtm_models=None,
              signalp_python=None, signalp_models=None):
    if Path(output).exists():
        raise FileExistsError(output)
    config = {'contract_version':'1', 'created_at':now(), 'tools':{}}
    for name, py, package, models in (
        ('IUPred2A',iupred_python,'iupred',None),
        ('DeepTMHMM2',dtm_python,'deeptmhmm2_predictor',dtm_models),
        ('SignalP',signalp_python,'signalp',signalp_models)):
        if not py:
            config['tools'][name] = {'status':'not_available','reason':'Local software/model installation not provided'}
            continue
        try:
            if name != 'IUPred2A' and not models:
                raise ValueError('Explicit model directory required; no implicit model download during run')
            probe = inventory(py,package,models)
            if name == 'DeepTMHMM2':
                required = ['esm2_t33_650M_UR50D.pt','esm2_t33_650M_UR50D-contact-regression.pt']
                if not all((Path(models)/'torch/hub/checkpoints'/n).is_file() for n in required):
                    raise ValueError('ESM2 weights absent in <model_dir>/torch/hub/checkpoints; finish model setup first')
            executable = Path(py).parent / {'IUPred2A':'iupred','DeepTMHMM2':'dtm2','SignalP':'signalp6'}[name]
            if not executable.is_file():
                raise ValueError('Local executable absent')
            config['tools'][name] = dict(probe,status='configured_not_execution_validated',
                executable=str(executable.absolute()), executable_sha256=file_hash(executable))
        except (OSError,ValueError,subprocess.SubprocessError) as exc:
            config['tools'][name] = {'status':'error','reason':str(exc)}
    write_json(output,config)
    return config


def _check_install(tool):
    actual = inventory(tool['python'],tool['package'],tool.get('model_dir'))
    if any(actual[k] != tool[k] for k in ('version','code_data_sha256','models')):
        raise ValueError('Installed code/data/model changed; create new configuration')
    if file_hash(tool['executable']) != tool['executable_sha256']:
        raise ValueError('Executable changed')


def run_predictors(reference_set, config_path, outdir, tools=TOOLS, resume=False,
                   timeout=3600, limit=None):
    from .harness import engine_hash, verify_bundle
    from .receiver_risk import read_reference
    setdir, root = Path(reference_set).resolve(),Path(outdir).resolve()
    definition = read_json(setdir/'set_definition.json')
    config = read_json(config_path)
    if root == setdir or setdir in root.parents:
        raise ValueError('Outputs must not modify the reference set')
    if not tools or len(set(tools)) != len(tools) or set(tools) - set(TOOLS):
        raise ValueError('Explicit supported tools, no duplicates')
    if timeout <= 0 or (limit is not None and limit <= 0):
        raise ValueError('Positive timeout/limit required')
    inputs = {'set_sha256':file_hash(setdir/'set_definition.json'),'reference_set_path':str(setdir),
        'config_sha256':file_hash(config_path),'tools':list(tools),'timeout':timeout,
        'limit':limit,'engine_sha256':engine_hash()}
    dest = root/digest(inputs)[:20]
    if dest.exists():
        if resume and verify_predictions(dest):
            for name in tools:
                if config['tools'].get(name,{}).get('status') == 'configured_not_execution_validated':
                    _check_install(config['tools'][name])
            return {'run_dir':str(dest),'summary':read_json(dest/'predictor_summary.json'),'execution':'reused'}
        raise FileExistsError('Preserve prior attempt; select a new output root')
    dest.mkdir(parents=True)
    rows, proteins = [], {}
    entries = list({e['accession']:e for e in definition['entries'] if e.get('protein_file')}.values())
    selected = entries[:limit] if limit else entries
    def load(entry):
        try:
            path = (setdir/entry['protein_file']).resolve()
            if not path.is_relative_to(setdir):
                raise ValueError('Reference path/hash conflict')
            p = read_reference(path,entry['protein_sha256'])
            p['sequence'] = sequence(p['sequence'])
            if p['accession'] != entry['accession']:
                raise ValueError('Accession conflict')
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,80}',p['accession']):
                raise ValueError('Unsafe accession identifier')
            return p,None
        except (OSError,ValueError,KeyError,TypeError) as exc:
            return None,{'accession':entry['accession'],'reason':'reference_validation: '+str(exc)}
    with ThreadPoolExecutor(max_workers=4) as pool:
        for p,error in pool.map(load,selected):
            if p:
                proteins[p['accession']] = p
            else:
                rows.extend(dict(error,software=name,status='error') for name in tools)
    records = []
    for name in tools:
        tool = config['tools'].get(name,{'status':'not_available','reason':'Tool not configured'})
        error = None
        if tool['status'] != 'configured_not_execution_validated':
            error = tool.get('reason','Tool unavailable')
        else:
            try:
                _check_install(tool)
            except (OSError,ValueError,subprocess.SubprocessError) as exc:
                error = str(exc)
        if error:
            rows.extend({'accession':acc,'software':name,'status':'not_available' if tool['status']=='not_available' else 'error',
                         'reason':error} for acc in proteins)
            continue
        folder = dest/name
        folder.mkdir()
        jobs = []
        for acc,p in proteins.items():
            sub = folder/acc
            sub.mkdir()
            fasta = sub/'input.fasta'
            fasta.write_text('>'+acc+'\n'+p['sequence']+'\n')
            jobs.append({'accession':acc,'sequence':p['sequence'],'fasta':str(fasta),
                         'output':str(sub/'scores.tsv' if name=='IUPred2A' else sub/'vendor')})
        parameters = {'execution':'local_only','device':'cpu','input_scope':'full_reference'}
        outcome = {}
        if name == 'IUPred2A':
            parameters.update(method='iupred2',mode='long',package_version=tool['version'])
            jobsfile = folder/'jobs.json'
            write_json(jobsfile,jobs)
            worker = Path(__file__).resolve().parents[1]/'local_iupred_worker.py'
            command = [tool['python'],str(worker),'--jobs',str(jobsfile),'--mode','long']
            # Scores already written before a timeout remain auditable; never fabricate the rest.
            with (folder/'stdout.log').open('w') as stdout, (folder/'stderr.log').open('w') as stderr:
                try:
                    proc = subprocess.run(command,stdout=stdout,stderr=stderr,timeout=timeout)
                    execution_error = None if proc.returncode == 0 else 'Worker exit '+str(proc.returncode)
                except subprocess.TimeoutExpired:
                    execution_error = 'Worker timeout'
            for line in (folder/'stdout.log').read_text().splitlines():
                try:
                    item = json.loads(line)
                    outcome[item['accession']] = item
                except (ValueError,KeyError):
                    continue
            write_json(folder/'invocation.json',{'command':command,'parameters':parameters,
                'software':tool,'completed_at':now(),'execution_error':execution_error})
        for job in jobs:
            acc, sub = job['accession'],Path(job['fasta']).parent
            try:
                if name == 'IUPred2A':
                    if outcome.get(acc,{}).get('status') != 'computed':
                        raise ValueError(outcome.get(acc,{}).get('reason',execution_error or 'No worker completion'))
                    raw = Path(job['output'])
                else:
                    if name == 'DeepTMHMM2':
                        parameters.update(simplify_io=True,batch_size=1)
                        command = [tool['executable'],job['fasta'],job['output'],'--model-dir',tool['model_dir'],
                                   '--device','cpu','--simplify-io','--batch-size','1']
                        raw = Path(job['output'])/'predicted_topologies.3line'
                    else:
                        parameters.update(organism='eukarya',mode='fast',format='none')
                        command = [tool['executable'],'--fastafile',job['fasta'],'--output_dir',job['output'],
                            '--model_dir',tool['model_dir'],'--organism','eukarya','--mode','fast','--format','none',
                            '--write_procs','1','--torch_num_threads','2']
                        raw = Path(job['output'])/'prediction_results.txt'
                    env = dict(os.environ,OMP_NUM_THREADS='2',MKL_NUM_THREADS='2')
                    if name == 'DeepTMHMM2':
                        env['TORCH_HOME'] = str(Path(tool['model_dir'])/'torch')
                    with (sub/'stdout.log').open('w') as stdout,(sub/'stderr.log').open('w') as stderr:
                        try:
                            proc = subprocess.run(command,stdout=stdout,stderr=stderr,timeout=timeout,env=env)
                            fail = None if proc.returncode==0 else 'Predictor exit '+str(proc.returncode)
                        except subprocess.TimeoutExpired:
                            fail = 'Predictor timeout'
                    write_json(sub/'invocation.json',{'command':command,'parameters':parameters,
                        'software':tool,'completed_at':now(),'execution_error':fail})
                    if fail:
                        raise ValueError(fail+'; see local stderr.log')
                rec = import_prediction(name,raw,job['fasta'],proteins[acc],tool['version'],dict(parameters))
                rec['execution_provenance'] = {'configuration_sha256':inputs['config_sha256'],
                    'code_data_sha256':tool['code_data_sha256'],'model_sha256':tool['models'],
                    'executable_sha256':tool['executable_sha256'],'executed_at':now(),
                    'version_semantics':'Installed package version; algorithm name recorded separately'}
                records.append(rec)
                rows.append({'accession':acc,'software':name,'status':'computed','reason':'',
                    'raw_output_sha256':rec['raw_output_sha256'],'version':tool['version']})
            except (OSError,ValueError,KeyError,TypeError) as exc:
                rows.append({'accession':acc,'software':name,'status':'error','reason':str(exc)})
    for name in tools:
        if config['tools'].get(name,{}).get('status')=='configured_not_execution_validated':
            _check_install(config['tools'][name])
    if file_hash(setdir/'set_definition.json') != inputs['set_sha256'] or file_hash(config_path)!=inputs['config_sha256']:
        raise ValueError('Input/config changed during execution')
    with ThreadPoolExecutor(max_workers=4) as pool:
        for p,error in pool.map(load,[e for e in selected if e['accession'] in proteins]):
            if error:
                raise ValueError('Reference changed during execution: '+error['reason'])
    complete = (definition['completeness']=='complete' and len(selected)==len(entries) and
                len(rows)==len(selected)*len(tools) and all(r['status']=='computed' for r in rows))
    summary = {'completeness':'complete' if complete else 'partial','requested_tools':list(tools),
        'set_reference_denominator':len(entries),'selected_reference_denominator':len(selected),
        'tool_reference_records':len(rows),'tool_states':{n:dict(Counter(r['status'] for r in rows if r['software']==n)) for n in tools},
        'unresolved_input_rows':definition['counts'].get('unresolved',0),'truncated':len(selected)<len(entries),
        'functional_accuracy':'not_estimated','claim_limit':'Predictor execution and annotation cross-checks are not SNIPR experimental accuracy.'}
    write_json(dest/'tool_evidence.json',records)
    write_json(dest/'predictor_summary.json',summary)
    tsv(dest/'predictor_status.tsv',rows,['accession','software','status','version','reason','raw_output_sha256'])
    write_json(dest/'inputs.json',inputs)
    files = {p.name:file_hash(p) for p in sorted(dest.iterdir()) if p.is_file()}
    nested = {str(p.relative_to(dest)):file_hash(p) for p in sorted(dest.rglob('*')) if p.is_file() and p.parent!=dest}
    write_json(dest/'manifest.json',{'state':'complete','artifacts':files,'nested_artifacts':nested,'inputs':inputs,'created_at':now()})
    return {'run_dir':str(dest),'summary':summary,'execution':'computed'}


def validate_local_tools(run_dir, reference_set, config, outdir, tools=TOOLS,
                         timeout=3600, interpro_cache=None, resume=False):
    """Predict reference set -> select explicit candidate references -> frozen overlay."""
    from .harness import verify_bundle
    from .receiver_risk import audit_risk
    source, root = Path(run_dir).resolve(),Path(outdir).resolve()
    if not verify_bundle(source):
        raise ValueError('Frozen screening source failed checksum verification')
    if root==source or source in root.parents:
        raise ValueError('Never overwrite frozen screening')
    predicted = run_predictors(reference_set,config,root/'prediction',tools,resume,timeout)
    pred_dir = Path(predicted['run_dir'])
    if not verify_predictions(pred_dir):
        raise ValueError('Predictor output verification failed')
    candidates = read_json(source/'candidates.json')
    accessions = {c.get('accession') or c.get('protein_id') for c in candidates}
    records = read_json(pred_dir/'tool_evidence.json')
    selected = [r for r in records if r['accession'] in accessions]
    path = root/('candidate-tool-evidence-'+digest(selected)[:20]+'.json')
    if path.exists():
        if not resume or read_json(path)!=selected:
            raise FileExistsError(path)
    else:
        write_json(path,selected)
    overlay = audit_risk(source,reference_set,root/'overlay',path,resume,interpro_cache)
    complete = predicted['summary']['completeness']=='complete' and overlay['summary']['completeness']=='complete'
    summary = {'completeness':'complete' if complete else 'partial',
        'prediction':predicted,'overlay':overlay,'records_for_candidate_references':len(selected),
        'records_without_candidate_reference':len(records)-len(selected),
        'source_manifest_sha256':file_hash(source/'manifest.json'),
        'functional_accuracy':'not_estimated','experimental_labels_used':False,
        'validation_meaning':'Execution, exact-sequence provenance and annotation cross-check; no held-out SNIPR functional benchmark.'}
    target = root/('validation-summary-'+digest(summary)[:20]+'.json')
    if not target.exists():
        write_json(target,summary)
    return summary
