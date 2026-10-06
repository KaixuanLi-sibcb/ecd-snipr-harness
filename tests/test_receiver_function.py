import copy
import csv
import tempfile
import unittest
from pathlib import Path

from test_design import project, protein, propose
from test_assays import observation
from test_screening import make_protein, make_set
from ecd_snipr.common import file_hash, read_json
from ecd_snipr.harness import run, verify_bundle
from ecd_snipr.observations import normalize_observation
from ecd_snipr.receiver_function import assess_receiver, audit_run, UNKNOWN
from ecd_snipr.screening import run_screening


def example():
    candidate = propose(protein())[0][0]
    candidate.update(fusion_sequence="MAAA" + candidate['sequence'] + "LLLRRR", scaffold_id="SYNTHETIC", scaffold_version="1")
    construct = {k: candidate[k] for k in ("protein_id", "start", "end", "fusion_sequence", "scaffold_id", "scaffold_version")}
    construct.update(construct_id="C1", antigen_sequence=candidate['sequence'], audit_status="sequence_consistent_not_functionally_validated")
    return candidate, {"C1": construct}


class ReceiverFunctionTests(unittest.TestCase):
    def test_standard_design_never_means_low_background(self):
        c = propose(protein())[0][0]
        c['screening_recommendation'] = 'standard_candidate'
        a = assess_receiver(c)
        self.assertEqual(a['functional_risk'], UNKNOWN)
        self.assertEqual(a['endpoints']['basal_activation']['evidence_state'], 'not_measured')
        self.assertIsNone(a['success_probability'])

    def test_no_risks_does_not_mean_safe(self):
        c = propose(protein())[0][0]
        c['risks'] = []
        a = assess_receiver(c)
        self.assertEqual(a['functional_risk'], UNKNOWN)
        self.assertTrue(all(x['state'] == 'not_established' for x in a['mechanism_review']))

    def test_generic_warnings_are_not_mechanistic_hits(self):
        c = propose(protein())[0][0]
        c['risks'] = [{'code': 'contextual_risk_evidence_missing'}, {'code': 'potential_n_glycosylation'}]
        a = assess_receiver(c)
        self.assertTrue(all(not x['signals'] for x in a['mechanism_review']))
        self.assertEqual(len(a['generic_warnings_not_failure_predictions']), 2)

    def test_name_blindness_and_no_mutation(self):
        c, _ = example()
        original = copy.deepcopy(c)
        a = assess_receiver(c)
        self.assertEqual(c, original)
        c.update(gene='UNSEEN-SYNTHETIC-GENE', accession='SYNTHETIC-OTHER', protein_id='OTHER')
        c['screening_recommendation'] = 'conditional_candidate'
        self.assertEqual(assess_receiver(c), a)

    def test_geometry_signal_not_functional_failure(self):
        c, _ = example()
        c['risks'] = [{'code': 'type_ii_attachment_orientation_change', 'detail': 'synthetic'}]
        a = assess_receiver(c)
        axis = next(x for x in a['mechanism_review'] if x['axis'] == 'attachment_geometry')
        self.assertEqual(axis['state'], 'review_signal_present')
        self.assertEqual(axis['functional_prediction'], 'not_made')
        self.assertEqual(a['endpoints']['basal_activation']['functional_risk'], UNKNOWN)

    def test_context_missing_outside_or_unresolved_not_positive_signal(self):
        for status, interpretation in [('missing', 'not_assessed'), ('unresolved', 'source_missing'),
                                       ('reported', 'reported_region_not_retained_not_proof_of_no_effect'),
                                       ('not_observed_in_context', 'only_this_context')]:
            c, _ = example()
            c['risks'] = []
            c['criteria_review'] = [{'kind': 'native_shedding', 'records': [{'effective_status': status, 'interpretation': interpretation}]}]
            axis = next(x for x in assess_receiver(c)['mechanism_review'] if x['axis'] == 'processing_junctions')
            self.assertFalse(axis['signals'])
            self.assertEqual(len(axis['context_records']), 1)

    def test_reported_native_shedding_is_review_not_basal_prediction(self):
        c, _ = example()
        c['criteria_review'] = [{'kind': 'native_shedding', 'records': [{'effective_status': 'reported', 'interpretation': 'requires_context_review', 'record': {'scope': 'native_protein', 'evidence': {'source': 'synthetic://source'}}}]}]
        a = assess_receiver(c)
        axis = next(x for x in a['mechanism_review'] if x['axis'] == 'processing_junctions')
        self.assertEqual(len(axis['signals']), 1)
        self.assertEqual(a['functional_risk'], UNKNOWN)

    def test_eligible_measurement_not_automatic_success(self):
        c, constructs = example()
        o = observation()
        o.update(endpoint='basal_activation', raw_value=0)
        a = assess_receiver(c, [normalize_observation(o, constructs)], constructs)
        e = a['endpoints']['basal_activation']
        self.assertEqual(e['evidence_state'], 'measured_in_recorded_context')
        self.assertEqual(e['measurement_count'], 1)
        self.assertEqual(e['functional_risk'], UNKNOWN)
        self.assertEqual(a['endpoints']['induced_response']['evidence_state'], 'not_measured')

    def test_backbone_difference_never_transfers_label(self):
        c, constructs = example()
        constructs['C1']['scaffold_version'] = '2'
        a = assess_receiver(c, [normalize_observation(observation(), constructs)], constructs)
        e = a['endpoints']['surface_expression']
        self.assertEqual(e['evidence_state'], 'related_or_unresolved_only')
        self.assertEqual(e['measurement_count'], 0)
        self.assertEqual(e['related_or_unresolved_records'][0]['link_scope'], 'antigen_only_not_equivalent_receptor')

    def test_fragment_or_fusion_difference_never_transfers_label(self):
        for field in ('antigen_sequence', 'fusion_sequence'):
            c, constructs = example()
            constructs['C1'][field] += 'A'
            a = assess_receiver(c, [normalize_observation(observation(), constructs)], constructs)
            self.assertEqual(a['endpoints']['surface_expression']['measurement_count'], 0)

    def test_unconfirmed_construct_audit_does_not_count(self):
        c, constructs = example()
        constructs['C1']['audit_status'] = 'requires_review'
        a = assess_receiver(c, [normalize_observation(observation(), constructs)], constructs)
        self.assertEqual(a['endpoints']['surface_expression']['measurement_count'], 0)

    def test_missing_and_y_n_do_not_become_failure_labels(self):
        for value in ('', 'Y', 'N', None):
            c, constructs = example()
            o = dict(observation(), raw_value=value)
            e = assess_receiver(c, [normalize_observation(o, constructs)], constructs)['endpoints']['surface_expression']
            self.assertEqual(e['evidence_state'], 'related_or_unresolved_only')
            self.assertEqual(e['functional_risk'], UNKNOWN)

    def test_conflicting_same_context_preserved(self):
        c, constructs = example()
        x, y = observation(), dict(observation(), observation_id='OBS2', raw_value=20)
        records = [normalize_observation(o, constructs) for o in (x, y)]
        e = assess_receiver(c, records, constructs)['endpoints']['surface_expression']
        self.assertEqual(e['evidence_state'], 'conflicting_measurements')
        self.assertEqual(e['measurement_count'], 2)
        self.assertEqual(len(e['exact_measurements']), 2)

    def test_batch_replicate_and_gate_not_pooled(self):
        for field in ('batch_id', 'replicate_id', 'condition_id'):
            c, constructs = example()
            x, y = observation(), dict(observation(), observation_id='OBS2', raw_value=20)
            y['context'][field] = 'OTHER'
            e = assess_receiver(c, [normalize_observation(o, constructs) for o in (x,y)], constructs)['endpoints']['surface_expression']
            self.assertEqual(e['evidence_state'], 'measured_in_recorded_context')
            self.assertEqual(len(e['measurement_contexts']), 2)
        y['gate_path'] = 'OTHER-GATE'
        e = assess_receiver(c, [normalize_observation(o, constructs) for o in (x,y)], constructs)['endpoints']['surface_expression']
        self.assertEqual(len(e['measurement_contexts']), 2)

    def test_downstream_screening_not_receiver_endpoint(self):
        c, constructs = example()
        o = dict(observation(), endpoint='antibody_screening_failure')
        a = assess_receiver(c, [normalize_observation(o, constructs)], constructs)
        self.assertTrue(all(e['measurement_count'] == 0 for e in a['endpoints'].values()))

    def test_synthetic_observation_never_biological_validation(self):
        c, constructs = example()
        c['synthetic_only'] = True
        e = assess_receiver(c, [normalize_observation(observation(), constructs)], constructs)['endpoints']['surface_expression']
        self.assertIn('synthetic_fixture_not_biological_validation', e['missing_information'])


class ReceiverIntegrationTests(unittest.TestCase):
    def test_generated_function_evidence_excluded_from_packages(self):
        from manage_skill import forbidden
        for name in ('receiver_function.tsv', 'receiver_mechanism_review.tsv', 'receiver_function_summary.json', 'RECEIVER_FUNCTION.md', 'source_link.json'):
            self.assertTrue(forbidden(Path('examples') / name))

    def test_screen_outputs_and_no_design_downgrade(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_set(root/'set', [make_protein('P98001')])
            result = run_screening(root/'set', root/'out', pilot=True)
            dest = Path(result['run_dir'])
            self.assertTrue(verify_bundle(dest))
            record = read_json(dest/'screening.json')[0]
            self.assertEqual(record['screening_recommendation'], 'standard_candidate')
            self.assertEqual(record['receiver_functional_risk'], UNKNOWN)
            with (dest/'receiver_function.tsv').open() as f:
                rows = list(csv.DictReader(f, delimiter='\t'))
            self.assertEqual(len(rows), 4)
            self.assertTrue(all(r['functional_risk'] == UNKNOWN for r in rows))
            self.assertIn('片段可设计不等于受体功能适配', (dest/'PI_SUMMARY.md').read_text())

    def test_overlay_preserves_design_and_input_hashes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = Path(run(project(), root/'source')['run_dir'])
            hashes = {p.name: file_hash(p) for p in original.iterdir() if p.is_file()}
            result = audit_run(original, root/'overlay')
            dest = Path(result['run_dir'])
            self.assertTrue(verify_bundle(dest))
            self.assertEqual(hashes, {p.name: file_hash(p) for p in original.iterdir() if p.is_file()})
            before, after = read_json(original/'candidates.json'), read_json(dest/'candidates.json')
            for a, b in zip(before, after):
                a.pop('receiver_function', None); b.pop('receiver_function', None)
                self.assertEqual(a, b)
            self.assertEqual(audit_run(original, root/'overlay', True)['execution'], 'verified_cache_hit')

    def test_overlay_refuses_original_nested_or_corrupt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = Path(run(project(), root/'source')['run_dir'])
            for dest in (original, original/'subdir'):
                with self.assertRaises(ValueError): audit_run(original, dest)
            result = audit_run(original, root/'overlay')
            (Path(result['run_dir'])/'receiver_function.tsv').write_text('broken')
            with self.assertRaises(FileExistsError): audit_run(original, root/'overlay', True)
            (original/'candidates.json').write_text('broken')
            with self.assertRaises(ValueError): audit_run(original, root/'new')


if __name__ == '__main__':
    unittest.main()
