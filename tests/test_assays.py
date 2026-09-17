import copy
import unittest
from test_design import protein
from ecd_snipr.observations import normalize_observation
from ecd_snipr.provenance import claims_for


def observation():
    return {"observation_id": "OBS1", "construct_id": "C1", "endpoint": "surface_expression", "raw_value": "80", "unit": "percent", "denominator": "live singlet receiver cells", "gate_path": "live/singlet/receiver", "channel": "synthetic-tag", "source_ref": {"workbook_sha256": "synthetic", "sheet": "test", "cell": "C2"},
            "context": {"batch_id": "B1", "host_cell": "synthetic", "sender_id": "none", "antibody_id": "not_applicable", "stimulus": "none", "condition_id": "COND1", "timepoint": "T1", "replicate_id": "R1"}}


CONSTRUCTS = {"C1": {"audit_status": "sequence_consistent_not_functionally_validated"}}


class AssayTests(unittest.TestCase):
    def test_percent_explicit(self):
        o = normalize_observation(observation(), CONSTRUCTS)
        self.assertEqual(o["normalized_value"], 0.8)
        self.assertEqual(o["eligibility"], "eligible_measurement_not_success_label")

    def test_blank_is_missing(self):
        o = observation()
        o["raw_value"] = ""
        result = normalize_observation(o, CONSTRUCTS)
        self.assertEqual(result["eligibility"], "missing")
        self.assertIsNone(result["normalized_value"])

    def test_Y_N_never_success(self):
        for value in ("Y", "N", "reconstruction", True):
            o = observation()
            o["raw_value"] = value
            self.assertEqual(normalize_observation(o, CONSTRUCTS)["eligibility"], "unresolved")

    def test_unknown_denominator(self):
        o = observation()
        del o["denominator"]
        self.assertEqual(normalize_observation(o, CONSTRUCTS)["eligibility"], "unresolved")

    def test_unit_not_inferred_from_value(self):
        o = observation()
        del o["unit"]
        self.assertIsNone(normalize_observation(o, CONSTRUCTS)["normalized_value"])

    def test_MFI_statistic_required(self):
        o = observation()
        o["unit"] = "fluorescence_au"
        self.assertIn("MFI_statistic_unconfirmed", normalize_observation(o, CONSTRUCTS)["issues"])

    def test_stimulus_required(self):
        o = observation()
        o["endpoint"] = "induced_response"
        self.assertIn("stimulated_condition_not_confirmed", normalize_observation(o, CONSTRUCTS)["issues"])

    def test_self_activation_requires_control(self):
        o = observation()
        o["endpoint"] = "basal_activation"
        o["context"]["stimulus"] = "unknown"
        self.assertIn("basal_condition_not_confirmed", normalize_observation(o, CONSTRUCTS)["issues"])

    def test_numeric_zero_is_not_missing(self):
        o = observation()
        o["raw_value"] = 0
        self.assertEqual(normalize_observation(o, CONSTRUCTS)["normalized_value"], 0)

    def test_nan_rejected(self):
        o = observation()
        o["raw_value"] = "nan"
        self.assertIsNone(normalize_observation(o, CONSTRUCTS)["normalized_value"])

    def test_different_batches_not_conflict(self):
        a, b = observation(), observation()
        b["observation_id"], b["raw_value"], b["context"]["batch_id"] = "OBS2", "20", "B2"
        records = [normalize_observation(x, CONSTRUCTS) for x in (a, b)]
        self.assertEqual(claims_for({"proteins": []}, [], [], records)[2], [])

    def test_same_context_conflict_not_overwritten(self):
        a, b = observation(), observation()
        b["observation_id"], b["raw_value"] = "OBS2", "20"
        records = [normalize_observation(x, CONSTRUCTS) for x in (a, b)]
        self.assertEqual(len(claims_for({"proteins": []}, [], [], records)[2]), 1)

    def test_unresolved_construct_not_training_label(self):
        self.assertEqual(normalize_observation(observation(), {})["eligibility"], "unresolved")


if __name__ == "__main__":
    unittest.main()
