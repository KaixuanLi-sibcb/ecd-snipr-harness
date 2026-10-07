import importlib.util
import contextlib
import io
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from test_design import ROOT
from test_screening import make_set, make_protein
from ecd_snipr.common import file_hash, read_json, write_json
from ecd_snipr.harness import verify_bundle
from ecd_snipr.local_predictors import verify_predictions
from ecd_snipr.predictor_import import import_prediction
from ecd_snipr.screening import run_screening
from manage_skill import forbidden


def helper(name):
    spec = importlib.util.spec_from_file_location('hpc_'+name, ROOT/'scripts/hpc'/f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


prepare = helper('prepare_jobs').prepare
reconcile = helper('reconcile')
finalize = helper('finalize').finalize
compat = helper('dtm2_api_compat')
signalp = helper('run_signalp')


class HPCWorkflowTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.setdir = self.root/'set'
        make_set(self.setdir,[make_protein('SYNTHETIC_A'),make_protein('SYNTHETIC_B')])
        self.jobs = self.root/'jobs.json'

    def simulate_results(self):
        jobs = read_json(self.jobs)
        for job in jobs['jobs']:
            folder = self.root/'raw'/'DeepTMHMM2'/job['accession']
            folder.mkdir(parents=True)
            protein = read_json(job['protein_path'])
            fasta = folder/'input.fasta'
            raw = folder/'predicted_topologies.3line'
            fasta.write_text('>'+job['accession']+'\n'+job['sequence']+'\n')
            raw.write_text('>'+job['accession']+'\n'+job['sequence']+'\n'+'O'*len(job['sequence'])+'\n')
            record = import_prediction('DeepTMHMM2',raw,fasta,protein,'synthetic-fixture',{'simplify_io':True})
            write_json(folder/'evidence.json',record)
            write_json(folder/'status.json',{'accession':job['accession'],'status':'computed',
                       'contract':{'jobs_sha256':file_hash(self.jobs)}})

    def reconcile(self):
        args = ['reconcile','--jobs',str(self.jobs),'--run',str(self.root/'raw'),
                '--outdir',str(self.root/'accepted'),'--tools','DeepTMHMM2']
        with patch.object(sys,'argv',args), contextlib.redirect_stdout(io.StringIO()):
            return reconcile.main()

    def test_exact_synthetic_pipeline_and_frozen_overlay(self):
        prepare(self.setdir,self.jobs)
        self.simulate_results()
        self.assertEqual(self.reconcile(),0)
        accepted = self.root/'accepted'
        self.assertTrue(verify_predictions(accepted))
        screen = run_screening(self.setdir,self.root/'screen')
        source = Path(screen['run_dir'])
        before = file_hash(source/'manifest.json')
        result = finalize(accepted,source,self.root/'overlay')
        self.assertTrue(verify_bundle(result['overlay']['run_dir']))
        self.assertEqual(before,file_hash(source/'manifest.json'))
        self.assertEqual(result['records_for_candidate_references'],2)
        self.assertEqual(result['functional_accuracy'],'not_estimated')

    def test_no_scaffold_needed_for_overlay(self):
        prepare(self.setdir,self.jobs)
        self.simulate_results()
        self.reconcile()
        screen = run_screening(self.setdir,self.root/'screen')
        self.assertTrue(finalize(self.root/'accepted',screen['run_dir'],self.root/'overlay'))

    def test_raw_tamper_rejected(self):
        prepare(self.setdir,self.jobs)
        self.simulate_results()
        self.reconcile()
        raw = self.root/'raw/DeepTMHMM2/SYNTHETIC_A/predicted_topologies.3line'
        raw.write_text('tampered\n')
        self.assertFalse(verify_predictions(self.root/'accepted'))
        with self.assertRaises(ValueError):
            finalize(self.root/'accepted')

    def test_resume_wrong_jobs_not_accepted(self):
        prepare(self.setdir,self.jobs)
        self.simulate_results()
        status = self.root/'raw/DeepTMHMM2/SYNTHETIC_A/status.json'
        record = read_json(status)
        record['contract']['jobs_sha256']='0'*64
        write_json(status,record)
        self.assertEqual(self.reconcile(),1)
        self.assertEqual(read_json(self.root/'accepted/summary.json')['states']['DeepTMHMM2']['integrity_error'],1)

    def test_one_missing_result_keeps_partial(self):
        prepare(self.setdir,self.jobs)
        self.simulate_results()
        (self.root/'raw/DeepTMHMM2/SYNTHETIC_B/status.json').unlink()
        self.assertEqual(self.reconcile(),1)
        summary=read_json(self.root/'accepted/summary.json')
        self.assertEqual(summary['completeness'],'partial')
        self.assertEqual(summary['states']['DeepTMHMM2']['not_run'],1)

    def test_unsupported_u_is_not_substituted(self):
        p=make_protein('SYNTHETIC_U')
        p['sequence']='U'+p['sequence'][1:]
        make_set(self.root/'u_set',[p])
        manifest=prepare(self.root/'u_set',self.jobs)
        self.assertEqual(manifest['eligible_total'],0)
        self.assertEqual(manifest['exclusions'][0]['status'],'unsupported_sequence')
        self.assertEqual(manifest['exclusions'][0]['sequence'],p['sequence'])

    def test_invalid_one_does_not_block_other(self):
        definition=read_json(self.setdir/'set_definition.json')
        (self.setdir/definition['entries'][0]['protein_file']).write_text('{}')
        manifest=prepare(self.setdir,self.jobs)
        self.assertEqual(manifest['eligible_total'],1)
        self.assertEqual(manifest['exclusions'][0]['status'],'reference_error')

    def test_duplicate_rows_do_not_duplicate_jobs(self):
        d=read_json(self.setdir/'set_definition.json')
        d['entries'].append(dict(d['entries'][0],duplicate_of='entry-00001'))
        write_json(self.setdir/'set_definition.json',d)
        manifest=prepare(self.setdir,self.jobs)
        self.assertEqual(manifest['source_input_rows'],3)
        self.assertEqual(manifest['scope_total'],2)

    def test_pilot_complete_eligible_still_partial_scope(self):
        prepare(self.setdir,self.jobs,pilot_limit=1)
        self.simulate_results()
        self.assertEqual(self.reconcile(),0)
        summary=read_json(self.root/'accepted/summary.json')
        self.assertTrue(summary['execution_complete_for_eligible'])
        self.assertEqual(summary['completeness'],'partial')
        self.assertEqual(summary['source_reference_denominator'],2)

    def test_partial_source_never_complete(self):
        d=read_json(self.setdir/'set_definition.json')
        d['completeness']='partial'
        write_json(self.setdir/'set_definition.json',d)
        prepare(self.setdir,self.jobs)
        self.simulate_results()
        self.reconcile()
        self.assertEqual(read_json(self.root/'accepted/summary.json')['completeness'],'partial')

    def test_unresolved_rows_preserved(self):
        d=read_json(self.setdir/'set_definition.json')
        d['entries'].append({'accession':None,'status':'ambiguous'})
        write_json(self.setdir/'set_definition.json',d)
        manifest=prepare(self.setdir,self.jobs)
        self.assertEqual(len(manifest['unresolved_source_rows']),1)
        self.simulate_results()
        self.reconcile()
        self.assertEqual(read_json(self.root/'accepted/summary.json')['completeness'],'partial')

    def test_preparation_cannot_overwrite_or_modify_source(self):
        prepare(self.setdir,self.jobs)
        with self.assertRaises(ValueError):
            prepare(self.setdir,self.jobs)
        with self.assertRaises(ValueError):
            prepare(self.setdir,self.setdir/'jobs.json')

    def test_adapter_single_batch_unpacking(self):
        self.assertEqual(compat.unpack_single_batch(['aaMMbb'],[[0.1,0.9]],[],6),
                         ('aaMMbb',[0.1,0.9],None))
        with self.assertRaises(ValueError):
            compat.unpack_single_batch(['aaM'],[[0.1]],[],6)
        with self.assertRaises(ValueError):
            compat.unpack_single_batch(['aa','bb'],[[0.1],[0.9]],[],2)

    def test_private_and_vendor_material_excluded(self):
        for path in ('models/a.pt','vendor/license.txt','envs/python','private_a/report.md',
                     'source.docx','software.tar.gz','prediction_jobs.json','pilot_jobs.json'):
            self.assertTrue(forbidden(Path(path)),path)

    def signalp_run(self, threads=1):
        config = self.root/'signalp-config.json'
        if not config.exists():
            write_json(config,{'tools':{'SignalP':{
                'executable':'/synthetic/signalp6','model_dir':'/synthetic/models',
                'version':'synthetic-fixture','code_data_sha256':'a'*64,
                'executable_sha256':'b'*64,'models':{}}}})
        args=['signalp','--jobs',str(self.jobs),'--config',str(config),
              '--outdir',str(self.root/'raw'),'--threads',str(threads)]
        def vendor(command,**kwargs):
            folder=Path(command[command.index('--output_dir')+1])
            folder.mkdir(parents=True)
            fasta=Path(command[command.index('--fastafile')+1]).read_text()
            ids=[line[1:] for line in fasta.splitlines() if line.startswith('>')]
            (folder/'prediction_results.txt').write_text('# ID\tPrediction\tCS Position\n'+
                ''.join(acc+'\tOTHER\t\n' for acc in ids))
            return subprocess.CompletedProcess(command,0)
        with patch.object(sys,'argv',args),patch.object(signalp,'_check_install'), \
             patch.object(signalp.subprocess,'run',side_effect=vendor) as runner, \
             contextlib.redirect_stdout(io.StringIO()):
            result=signalp.main()
        return result,runner.call_count

    def test_signalp_synthetic_batch_separation_and_resume(self):
        prepare(self.setdir,self.jobs)
        self.assertEqual(self.signalp_run(),(0,1))
        for job in read_json(self.jobs)['jobs']:
            record=read_json(self.root/'raw/SignalP'/job['accession']/'evidence.json')
            reconcile.verify_record(record,read_json(job['protein_path']))
        self.assertEqual(self.signalp_run(),(0,0))

    def test_signalp_changed_contract_preserves_status(self):
        prepare(self.setdir,self.jobs)
        self.signalp_run()
        state=self.root/'raw/SignalP/SYNTHETIC_A/status.json'
        checksum=file_hash(state)
        with self.assertRaises(ValueError):
            self.signalp_run(threads=2)
        self.assertEqual(checksum,file_hash(state))

    def test_signalp_changed_row_preserves_prior_result(self):
        prepare(self.setdir,self.jobs)
        self.signalp_run()
        folder=self.root/'raw/SignalP/SYNTHETIC_A'
        checksum=file_hash(folder/'status.json')
        (folder/'prediction_results.txt').write_text('tampered\n')
        with self.assertRaises(ValueError):
            self.signalp_run()
        self.assertEqual(checksum,file_hash(folder/'status.json'))

    def test_signalp_reimport_rejects_wrong_batch_row(self):
        prepare(self.setdir,self.jobs)
        self.signalp_run()
        folder=self.root/'raw/SignalP/SYNTHETIC_A'
        record=read_json(folder/'evidence.json')
        (folder/'prediction_results.txt').write_text('# ID\tPrediction\tCS Position\nSYNTHETIC_B\tOTHER\t\n')
        with self.assertRaises(ValueError):
            reconcile.verify_record(record,read_json(read_json(self.jobs)['jobs'][0]['protein_path']))

    def test_scheduler_example_refuses_unallocated_execution(self):
        import os
        env=dict(os.environ)
        env.pop('JOB_ID',None)
        result=subprocess.run(['bash',str(ROOT/'examples/hpc_worker.sh')],env=env,capture_output=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn(b'qsub',result.stderr)


if __name__=='__main__':
    unittest.main()
