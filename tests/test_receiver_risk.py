import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import subprocess

from test_design import feature
from test_screening import make_protein, make_set
from ecd_snipr.common import digest, file_hash, read_json, write_json
from ecd_snipr.design import propose
from ecd_snipr.harness import verify_bundle
from ecd_snipr.receiver_risk import assess_risk, audit_risk, run_sequence_tools, read_reference
from ecd_snipr.sequence_tools import descriptors, sequence_hash, validate_prediction
from ecd_snipr.screening import run_screening


def case():
    p = make_protein('SYNTHETIC')
    return p, propose(p)[0][0]


class RiskTests(unittest.TestCase):
    def test_no_specific_signal_is_not_safe(self):
        p,c = case()
        r = assess_risk(c,p)
        self.assertEqual(r['review_priority'],'no_specific_signal_detected')
        self.assertEqual(r['functional_risk'],'undetermined_not_low_risk')
        self.assertIsNone(r['success_probability'])

    def test_optional_context_missing_does_not_block_design(self):
        p,c = case()
        r = assess_risk(c,p)
        self.assertEqual(r['assessment_status'],'annotation_and_sequence_evaluated')
        self.assertEqual(c['design_status'],'needs_review')

    def test_name_and_lab_outcome_blind(self):
        p,c = case()
        p['subunit_comments']=[{'text':'Homodimer.'}]
        a = assess_risk(c,p)
        p['gene']=c['gene']='DIFFERENT_UNSEEN_GENE'
        p['accession']=c['accession']='ANOTHER_ACCESSION'
        p['self_activation']=c['self_activation']=True
        p['failure_reason']='known failure'
        self.assertEqual(assess_risk(c,p),a)

    def test_no_input_mutation(self):
        p,c = case()
        original=copy.deepcopy((p,c))
        assess_risk(c,p)
        self.assertEqual((p,c),original)

    def test_whole_protein_homodimer_only_context(self):
        p,c = case()
        p['subunit_comments']=[{'text':'Homodimer.'}]
        r=assess_risk(c,p)
        self.assertEqual(r['review_priority'],'context_dependent_review')
        self.assertIn('native_self_association_context',r['reason_codes'])

    def test_homophilic_native_full_ecd_elevated_not_measured(self):
        p,c=case()
        p['function_comments']=[{'text':'Mediates homophilic binding.'}]
        r=assess_risk(c,p)
        self.assertEqual(r['review_priority'],'elevated_review_priority')
        self.assertEqual(r['functional_confidence'],'not_established')

    def test_domain_fragment_does_not_inherit_whole_ecd_association(self):
        p,c=case()
        p['function_comments']=[{'text':'Mediates homophilic binding.'}]
        c['antigen_form_type']='domain_fragment'
        self.assertEqual(assess_risk(c,p)['review_priority'],'context_dependent_review')

    def test_self_ligand_annotation_alias_is_general_not_gene_rule(self):
        for term in ('Self-ligand receptor','Self ligand receptor'):
            p,c=case(); p['function_comments']=[{'text':term}]
            self.assertEqual(assess_risk(c,p)['review_priority'],'elevated_review_priority')

    def test_negated_self_ligand_annotation_is_unresolved(self):
        p,c=case(); p['function_comments']=[{'text':'Not a self-ligand receptor.'}]
        self.assertEqual(assess_risk(c,p)['review_priority'],'no_specific_signal_detected')

    def test_intracellular_association_not_transferred(self):
        p,c=case()
        p['subunit_comments']=[{'text':'Homodimer via the transmembrane domain.'}]
        r=assess_risk(c,p)
        self.assertEqual(r['review_priority'],'no_specific_signal_detected')
        self.assertTrue(r['excluded_or_unresolved_evidence'])

    def test_negation_not_positive(self):
        for field,text in [('subunit_comments','Does not form homodimers.'),('ptm_comments','Not shed by ADAM10.'),('function_comments','Does not bind sialic acid.')]:
            p,c=case(); p[field]=[{'text':text}]
            r=assess_risk(c,p)
            self.assertEqual(r['review_priority'],'no_specific_signal_detected')
            self.assertTrue(r['excluded_or_unresolved_evidence'])

    def test_retained_processing_site_not_native_text_only(self):
        p,c=case()
        p['features'].append(feature('site',30,30,name='Cleavage by a protease'))
        r=assess_risk(c,p)
        self.assertEqual(r['review_priority'],'elevated_review_priority')
        self.assertEqual(r['signals'][0]['evidence']['reference_interval'],[30,30])

    def test_processing_outside_fragment_not_counted(self):
        p,c=case()
        p['features'].append(feature('site',90,90,name='Cleavage by a protease'))
        self.assertEqual(assess_risk(c,p)['review_priority'],'no_specific_signal_detected')

    def test_unlocalized_shedding_only_context(self):
        p,c=case(); p['ptm_comments']=[{'text':'The ectodomain is shed by ADAM10.'}]
        self.assertEqual(assess_risk(c,p)['review_priority'],'context_dependent_review')

    def test_regulating_another_proteins_shedding_not_self_processing(self):
        p,c=case(); p['function_comments']=[{'text':'Regulates ectodomain shedding of another receptor.'}]
        self.assertEqual(assess_risk(c,p)['review_priority'],'no_specific_signal_detected')

    def test_adam_name_alone_not_shedding_evidence(self):
        p,c=case(); p['ptm_comments']=[{'text':'Interacts with ADAM10.'}]
        self.assertEqual(assess_risk(c,p)['review_priority'],'no_specific_signal_detected')

    def test_retained_interchain_disulfide_not_cysteine_count(self):
        p,c=case(); p['features'].append(feature('disulfide',30,30,name='Interchain'))
        r=assess_risk(c,p)
        self.assertIn('retained_interchain_disulfide',r['reason_codes'])
        self.assertEqual(r['review_priority'],'elevated_review_priority')

    def test_generic_glyco_and_odd_cysteine_do_not_escalate(self):
        p,c=case(); c['risks']=[{'code':'potential_n_glycosylation'},{'code':'free_thiol_odd_cysteine'}]
        self.assertEqual(assess_risk(c,p)['review_priority'],'no_specific_signal_detected')

    def test_source_and_sequence_mismatch_block_assessment(self):
        p,c=case(); p['sequence']='A'*len(p['sequence'])
        self.assertEqual(assess_risk(c,p)['assessment_status'],'reference_conflict')
        p,c=case(); p['evidence']={}
        self.assertEqual(assess_risk(c,p)['review_priority'],'not_assessed')

    def test_invalid_or_blocked_candidate_not_evaluated(self):
        p,c=case(); c['sequence']='INVALID-X'
        self.assertEqual(assess_risk(c,p)['assessment_status'],'invalid_sequence')
        p,c=case(); c['design_status']='blocked'
        self.assertEqual(assess_risk(c,p)['review_priority'],'not_assessed')


class SequenceTests(unittest.TestCase):
    def test_kd_and_entropy_known_extremes(self):
        d=descriptors('L'*40)
        self.assertAlmostEqual(d['gravy'],3.8)
        self.assertEqual(d['hydrophobic_segments'][0]['end'],40)
        self.assertEqual(d['low_complexity_segments'][0]['extreme_value'],0)
        self.assertFalse(d['stp_enriched_segments'])
        self.assertTrue(descriptors('STP'*20)['stp_enriched_segments'])

    def test_short_sequence_has_no_forced_window(self):
        d=descriptors('ACDE')
        self.assertEqual(d['hydrophobic_segments'],[])
        self.assertEqual(d['low_complexity_segments'],[])

    def test_sequence_alert_never_elevated(self):
        _,c=case(); c['sequence']='L'*40
        self.assertEqual(assess_risk(c)['review_priority'],'sequence_alert_only')

    def test_optional_biopython_explicit_state(self):
        d=descriptors('ACDEFGHIKLMNPQRSTVWY'*2,True)
        self.assertIn(d['biopython']['status'],{'computed','not_available'})
        if d['biopython']['status']=='computed':
            self.assertEqual(d['biopython']['kd_crosscheck'],'pass')

    def prediction(self):
        _,c=case()
        r={'software':'IUPred3','version':'test','candidate_id':c['candidate_id'],
           'sequence_sha256':sequence_hash(c['sequence']),'input_scope':'candidate_fragment',
           'coordinate_system':'1-based-inclusive','parameters':{},'source':'synthetic fixture',
           'raw_output_sha256':'a'*64,'scores':[0.8]*len(c['sequence'])}
        return c,r

    def test_disorder_not_functional_risk(self):
        c,r=self.prediction(); result=assess_risk(c,predictions=[r])
        self.assertEqual(result['review_priority'],'no_specific_signal_detected')
        self.assertEqual(result['disorder_descriptors'][0]['fraction_above_0_5'],1)

    def test_predictor_contract_rejects_hash_scope_nan_and_coordinates(self):
        c,r=self.prediction()
        for key,value in [('sequence_sha256','wrong'),('input_scope','whole_protein'),('coordinate_system','0-based'),('scores',[float('nan')]*len(c['sequence']))]:
            bad=dict(r,**{key:value})
            with self.assertRaises(ValueError): validate_prediction(bad,c)
        bad=dict(r,software='DeepTMHMM',regions=[{'kind':'transmembrane','start':0,'end':20}])
        with self.assertRaises(ValueError): validate_prediction(bad,c)

    def test_predicted_tm_requires_review_not_failure(self):
        c,r=self.prediction(); r.update(software='DeepTMHMM',regions=[{'kind':'transmembrane','start':20,'end':35}]); r.pop('scores')
        a=assess_risk(c,predictions=[r])
        self.assertEqual(a['review_priority'],'context_dependent_review')
        self.assertEqual(a['functional_confidence'],'not_established')

    def test_rejected_prediction_remains_visible(self):
        c,r=self.prediction(); r['sequence_sha256']='wrong'
        a=assess_risk(c,predictions=[r])
        self.assertFalse(a['tool_evidence'])
        self.assertEqual(a['excluded_or_unresolved_evidence'][0]['reason'],'rejected_tool_evidence')


class RiskIntegrationTests(unittest.TestCase):
    def test_blocked_file_read_is_bounded(self):
        with patch('ecd_snipr.receiver_risk.subprocess.run',side_effect=subprocess.TimeoutExpired('reader',5)):
            with self.assertRaisesRegex(ValueError,'reference_read_timeout'): read_reference('synthetic','a'*64)

    def test_reference_hash_is_checked_on_exact_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'p.json'; write_json(path,{'sequence':'ACDE'})
            self.assertEqual(read_reference(path,file_hash(path)),{'sequence':'ACDE'})
            with self.assertRaisesRegex(ValueError,'hash conflict'): read_reference(path,'a'*64)

    def setup_bundle(self,root):
        p=make_protein('P98010'); p['ptm_comments']=[{'text':'The ectodomain is shed by ADAM10.'}]
        make_set(root/'set',[p])
        return Path(run_screening(root/'set',root/'screen',pilot=True)['run_dir'])

    def test_screen_exports_risk_without_design_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            run=self.setup_bundle(Path(tmp))
            self.assertTrue(verify_bundle(run))
            c=read_json(run/'candidates.json')[0]
            self.assertEqual(c['screening_recommendation'],'standard_candidate')
            self.assertEqual(c['receiver_risk']['review_priority'],'context_dependent_review')
            self.assertTrue((run/'receiver_risk.tsv').is_file())

    def test_overlay_immutable_reusable_and_reference_hash_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); run=self.setup_bundle(root); before=file_hash(run/'manifest.json')
            a=audit_risk(run,root/'set',root/'risk')
            self.assertTrue(verify_bundle(a['run_dir']))
            self.assertEqual(audit_risk(run,root/'set',root/'risk',resume=True)['execution'],'verified_cache_hit')
            self.assertEqual(file_hash(run/'manifest.json'),before)
            path=root/'set/proteins/P98010.json'; path.write_text('{}')
            with self.assertRaises(FileExistsError): audit_risk(run,root/'set',root/'risk',resume=True)
            b=audit_risk(run,root/'set',root/'partial')
            self.assertEqual(b['summary']['completeness'],'partial')
            self.assertEqual(b['summary']['usable_candidate_denominator'],1)

    def test_overlay_cannot_write_into_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); run=self.setup_bundle(root)
            with self.assertRaises(ValueError): audit_risk(run,root/'set',run/'risk')

    def test_sequence_tools_verify_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); run=self.setup_bundle(root)
            result=run_sequence_tools(run,root/'tools')
            self.assertTrue(verify_bundle(result['run_dir']))
            with self.assertRaises(FileExistsError): run_sequence_tools(run,root/'tools')

    def test_generated_outputs_not_packaged(self):
        from manage_skill import forbidden
        for name in ('receiver_risk.tsv','receiver_risk_signals.tsv','receiver_risk_summary.json','risk_policy_snapshot.json','sequence_descriptors.json'):
            self.assertTrue(forbidden(Path('examples')/name))
