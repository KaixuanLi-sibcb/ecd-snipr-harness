import tempfile
import json
import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch

from test_receiver_risk import case
from test_screening import make_set, make_protein
from ecd_snipr.common import read_json, write_json, file_hash
from ecd_snipr.local_predictors import configure, run_predictors, verify_predictions, inventory
from ecd_snipr.predictor_import import import_prediction
from ecd_snipr.core_evidence import check_core


class PredictorFormatTests(unittest.TestCase):
    def test_signalp_official_distribution_name(self):
        with tempfile.TemporaryDirectory() as temp:
            module = Path(temp)/'signalp'/'__init__.py'
            module.parent.mkdir()
            module.write_text('# synthetic package\n')
            response = subprocess.CompletedProcess([],0,json.dumps({'version':'6.0+h',
                'distribution':'signalp6','path':str(module)}),'')
            with patch('ecd_snipr.local_predictors._run',return_value=response) as runner:
                result = inventory('/synthetic/python','signalp')
            self.assertIn("d='signalp6'",runner.call_args.args[0][2])
            self.assertEqual(result['distribution'],'signalp6')
            self.assertEqual(result['package'],'signalp')

    def run_import(self, software, text, parameters=None):
        p,c = case()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'input.fasta').write_text('>'+p['accession']+'\n'+p['sequence']+'\n')
            (root/'raw').write_text(text.replace('ACCESSION',p['accession']).replace('SEQUENCE',p['sequence']))
            return import_prediction(software,root/'raw',root/'input.fasta',p,'synthetic-fixture',parameters or {}),p,c

    def test_signalp_explicit_no_signal(self):
        r,p,c = self.run_import('SignalP','# ID\tPrediction\tOTHER\tSP\tCS Position\nACCESSION\tOTHER\t0.99\t0.01\t\n')
        self.assertEqual(r['regions'],[])
        self.assertEqual(check_core(c,p,predictions=[r])['coverage']['SignalP'],'evaluated')

    def test_signalp_no_signal_is_not_blank_negative(self):
        with self.assertRaises(ValueError):
            self.run_import('SignalP','##gff-version 3\n')

    def test_signalp_cleavage_adjacent(self):
        r,_,_ = self.run_import('SignalP','# ID\tPrediction\tCS Position\nACCESSION\tSP\tCS pos: 10-11.\n')
        self.assertEqual(r['regions'][0]['end'],10)
        with self.assertRaises(ValueError):
            self.run_import('SignalP','# ID\tPrediction\tCS Position\nACCESSION\tSP\tCS pos: 10-12.\n')

    def test_signalp_unknown_class_rejected(self):
        with self.assertRaises(ValueError):
            self.run_import('SignalP','# ID\tPrediction\tCS Position\nACCESSION\tUNKNOWN\t\n')

    def test_signalp_identity_rejected(self):
        with self.assertRaises(ValueError):
            self.run_import('SignalP','# ID\tPrediction\tCS Position\nWRONG\tOTHER\t\n')

    def test_dtm2_explicit_distinct_version(self):
        p,c=case()
        r,_,_ = self.run_import('DeepTMHMM2','>ACCESSION\nSEQUENCE\n'+'O'*len(p['sequence'])+'\n',{'simplify_io':True})
        self.assertEqual(r['software'],'DeepTMHMM2')
        result = check_core(c,p,predictions=[r])
        self.assertEqual(result['coverage']['DeepTMHMM'],'evaluated')
        self.assertEqual(result['coverage']['disorder'],'not_run')

    def test_dtm2_needs_explicit_simplified_mode(self):
        p,c=case()
        with self.assertRaises(ValueError):
            self.run_import('DeepTMHMM2','>ACCESSION\nSEQUENCE\n'+'O'*len(p['sequence'])+'\n')

    def test_dtm2_new_labels_not_silently_legacy(self):
        p,c=case()
        labels = 'R'*len(p['sequence'])
        r,_,_ = self.run_import('DeepTMHMM2','>ACCESSION\nSEQUENCE\n'+labels+'\n',{'simplify_io':True})
        self.assertEqual(r['regions'][0]['kind'],'reentrant')
        self.assertIn('predicted_exclusion_region_retained',check_core(c,p,predictions=[r])['reason_codes'])
        with self.assertRaises(ValueError):
            self.run_import('DeepTMHMM','>ACCESSION\nSEQUENCE\n'+labels+'\n')


class RunnerTests(unittest.TestCase):
    def test_worker_single_failure_and_derived_tamper(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            make_set(root/'set',[make_protein('SYNTHETIC_A'),make_protein('SYNTHETIC_B')])
            configure(root/'config.json')
            config=read_json(root/'config.json')
            config['tools']['IUPred2A']={'status':'configured_not_execution_validated','python':'synthetic',
                'version':'fixture','code_data_sha256':'a'*64,'models':{},'executable_sha256':'b'*64}
            write_json(root/'config.json',config)
            original_run=subprocess.run
            def worker(command,**kwargs):
                if '--jobs' not in command:
                    return original_run(command,**kwargs)
                jobs=read_json(command[command.index('--jobs')+1])
                good=jobs[0]
                Path(good['output']).write_text('\n'.join(f'{i}\t{aa}\t0.25' for i,aa in enumerate(good['sequence'],1)))
                kwargs['stdout'].write('{"accession":"SYNTHETIC_A","status":"computed"}\n')
                kwargs['stdout'].write('{"accession":"SYNTHETIC_B","status":"error","reason":"fixture error"}\n')
                return subprocess.CompletedProcess(command,0)
            with patch('ecd_snipr.local_predictors._check_install'),patch('ecd_snipr.local_predictors.subprocess.run',worker):
                r=run_predictors(root/'set',root/'config.json',root/'out',tools=['IUPred2A'])
            self.assertEqual(r['summary']['tool_states']['IUPred2A'],{'computed':1,'error':1})
            self.assertEqual(r['summary']['completeness'],'partial')
            dest=Path(r['run_dir']); self.assertTrue(verify_predictions(dest))
            records=read_json(dest/'tool_evidence.json'); records[0]['scores'][0]=0.9
            write_json(dest/'tool_evidence.json',records)
            manifest=read_json(dest/'manifest.json'); manifest['artifacts']['tool_evidence.json']=file_hash(dest/'tool_evidence.json')
            write_json(dest/'manifest.json',manifest)
            self.assertFalse(verify_predictions(dest))

    def test_timeout_does_not_become_negative(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); make_set(root/'set',[make_protein('SYNTHETIC_A')]); configure(root/'config.json')
            config=read_json(root/'config.json')
            config['tools']['IUPred2A']={'status':'configured_not_execution_validated','python':'synthetic',
                'version':'fixture','code_data_sha256':'a'*64,'models':{},'executable_sha256':'b'*64}
            write_json(root/'config.json',config)
            original_run=subprocess.run
            def worker(command,**kwargs):
                if '--jobs' in command:
                    raise subprocess.TimeoutExpired('fixture',1)
                return original_run(command,**kwargs)
            with patch('ecd_snipr.local_predictors._check_install'),patch('ecd_snipr.local_predictors.subprocess.run',worker):
                r=run_predictors(root/'set',root/'config.json',root/'out',tools=['IUPred2A'])
            self.assertEqual(r['summary']['tool_states']['IUPred2A'],{'error':1})
            self.assertEqual(read_json(Path(r['run_dir'])/'tool_evidence.json'),[])

    def test_missing_install_is_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            config=configure(root/'config.json')
            self.assertEqual(config['tools']['SignalP']['status'],'not_available')
            with self.assertRaises(FileExistsError):
                configure(root/'config.json')

    def test_unavailable_tool_partial_not_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); make_set(root/'set',[make_protein('SYNTHETIC_A')]); configure(root/'config.json')
            r=run_predictors(root/'set',root/'config.json',root/'out')
            self.assertEqual(r['summary']['completeness'],'partial')
            self.assertEqual(read_json(Path(r['run_dir'])/'tool_evidence.json'),[])
            self.assertTrue(verify_predictions(r['run_dir']))
            again=run_predictors(root/'set',root/'config.json',root/'out',resume=True)
            self.assertEqual(again['execution'],'reused')
            Path(r['run_dir'],'predictor_status.tsv').write_text('tampered')
            self.assertFalse(verify_predictions(r['run_dir']))
            with self.assertRaises(FileExistsError):
                run_predictors(root/'set',root/'config.json',root/'out',resume=True)

    def test_reference_error_is_recorded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); make_set(root/'set',[make_protein('SYNTHETIC_A')]); configure(root/'config.json')
            d=read_json(root/'set/set_definition.json')
            Path(root/'set',d['entries'][0]['protein_file']).write_text('{}')
            r=run_predictors(root/'set',root/'config.json',root/'out',tools=['IUPred2A'])
            self.assertEqual(r['summary']['tool_states']['IUPred2A'],{'error':1})

    def test_unsupported_tool_or_negative_timeout(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); make_set(root/'set',[make_protein('SYNTHETIC_A')]); configure(root/'config.json')
            for kwargs in ({'tools':['Wrong']},{'timeout':0},{'limit':0}):
                with self.assertRaises(ValueError):
                    run_predictors(root/'set',root/'config.json',root/'out',**kwargs)


if __name__ == '__main__':
    unittest.main()
