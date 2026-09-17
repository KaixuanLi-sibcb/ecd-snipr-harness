import tempfile
import unittest
from pathlib import Path

from test_design import codes, feature, project, protein, scaffold, approve
from ecd_snipr.contracts import validate_project
from ecd_snipr.design import propose
from ecd_snipr.harness import run, verify_bundle
from ecd_snipr.uniprot import normalize


def secreted():
    p = protein()
    p.update(topology="secreted", location="secreted", candidate_regions=[],
             features=[feature("signal_peptide", 1, 10), feature("propeptide", 11, 15), feature("chain", 16, 95)])
    return p


def consideration(kind="native_shedding", **extra):
    return dict(kind=kind, status="reported", scope="native_protein", reason="Synthetic report", evidence=protein()["evidence"], context={"system": "synthetic_native_context"}, **extra)


class CriteriaTests(unittest.TestCase):
    def test_default_full_form_only(self):
        p = protein()
        p.pop("candidate_policy")
        p["candidate_regions"] = []
        self.assertEqual(len(propose(p)[0]), 1)

    def test_noncanonical_not_discarded(self):
        p = protein()
        p["isoform"] = "SYNTHETIC-001-2"
        self.assertNotEqual(propose(p)[0][0]["design_status"], "blocked")

    def test_secreted_chain_not_surface_target(self):
        c = propose(secreted())[0][0]
        self.assertEqual(c["antigen_form_type"], "mature_secreted")
        self.assertEqual(c["start"], 16)
        self.assertIn("secreted_antigen_scope_extension", codes(c))
        self.assertNotEqual(c["design_status"], "blocked")
        self.assertEqual(c["functional_status"], "not_experimentally_validated")

    def test_no_chain_no_signal_to_end_guess(self):
        p = secreted()
        p["features"] = p["features"][:2]
        self.assertEqual(propose(p)[0], [])

    def test_multiple_chains_not_independent_successes(self):
        p = secreted()
        p["features"][-1] = feature("chain", 16, 40)
        p["features"].append(feature("chain", 50, 95))
        for c in propose(p)[0]:
            self.assertIn("multiple_processed_chains_dependency_unknown", codes(c))

    def test_secreted_tm_conflict(self):
        p = secreted()
        p["features"].append(feature("transmembrane", 70, 90))
        self.assertEqual(propose(p)[0][0]["design_status"], "blocked")

    def test_propeptide_not_retained(self):
        p = secreted()
        p["features"][-1]["start"] = 11
        self.assertIn("retained_propeptide", codes(propose(p)[0][0]))

    def test_omega_retained_not_minus_one(self):
        p = protein()
        p.update(topology="gpi", features=[feature("signal_peptide", 1, 10), feature("chain", 11, 76), feature("gpi_attachment_site", 76, 76), feature("gpi_signal", 77, 95)], candidate_regions=[])
        c = propose(p)[0][0]
        self.assertEqual((c["start"], c["end"]), (11, 76))
        self.assertEqual(c["sequence"], p["sequence"][10:76])
        self.assertNotEqual(c["design_status"], "blocked")

    def test_gpi_omega_unknown_blocks_mature_claim(self):
        p = protein()
        p["topology"] = "gpi"
        self.assertIn("mature_gpi_boundary_not_supported", codes(propose(p)[0][0]))

    def test_gpi_signal_retention_without_signal_feature(self):
        p = protein()
        p["topology"] = "gpi"
        p["features"].append(feature("gpi_attachment_site", 70, 70))
        self.assertIn("retained_gpi_c_terminal_signal", codes(propose(p)[0][0]))

    def test_shed_not_from_motif(self):
        p = protein()
        r = dict(feature("region", 11, 65), antigen_form_type="shed_product", rationale="Test", release_evidence=dict(p["evidence"], kind="prediction"), release_context="synthetic")
        p["candidate_regions"] = [r]
        self.assertIn("shed_release_evidence_missing", codes(propose(p)[0][-1]))
        r["release_evidence"] = p["evidence"]
        c = propose(p)[0][-1]
        self.assertNotEqual(c["design_status"], "blocked")
        self.assertIn("shed_product_tethering_unvalidated", codes(c))

    def test_same_bounds_distinct_form_not_silently_merged(self):
        p = protein()
        p["candidate_regions"] = [dict(feature("region", 11, 76), antigen_form_type="shed_product", rationale="Test", release_evidence=p["evidence"], release_context="synthetic")]
        matches = [c for c in propose(p)[0] if (c["start"], c["end"]) == (11, 76)]
        self.assertEqual(len(matches), 2)
        self.assertNotEqual(matches[0]["candidate_id"], matches[1]["candidate_id"])

    def test_unknown_not_intracellular_failure(self):
        p = protein()
        p["topology"] = "unknown"
        self.assertIn("topology_unknown", codes(propose(p)[0][0]))
        p["topology"] = "intracellular"
        cs, issues = propose(p)
        self.assertEqual(cs, [])
        self.assertEqual(issues[0]["code"], "outside_extracellular_antigen_scope")

    def test_blank_risk_not_low_risk(self):
        c = propose(protein())[0][0]
        self.assertEqual({r["status"] for r in c["criteria_review"]}, {"missing"})
        self.assertEqual(c["risk_review_status"], "incomplete")

    def test_native_shedding_not_high_basal_label(self):
        p = protein()
        p["considerations"] = [consideration()]
        c = propose(p)[0][0]
        self.assertIn("native_shedding_context_review", codes(c))
        self.assertEqual(c["functional_status"], "not_experimentally_validated")
        self.assertIn("not_a_prediction", c["criteria_review"][1]["records"][0]["interpretation"])

    def test_source_free_consideration_unresolved(self):
        p = protein()
        r = consideration()
        r["evidence"] = {}
        p["considerations"] = [r]
        self.assertEqual(propose(p)[0][0]["criteria_review"][1]["status"], "unresolved")

    def test_no_pooling_positive_negative_contexts(self):
        p = protein()
        a, b = consideration(), consideration()
        b.update(status="not_observed_in_context", context={"system": "other"})
        p["considerations"] = [a, b]
        self.assertEqual(propose(p)[0][0]["criteria_review"][1]["status"], "mixed_context_records")

    def test_gene_name_not_automatic_b_cell_failure(self):
        p = protein()
        p["gene"] = "CD40"
        c = propose(p)[0][0]
        self.assertEqual(c["criteria_review"][-1]["status"], "missing")
        p["considerations"] = [dict(consideration("culture_interference"), scope="assay_context", context={"construct_id": "test", "condition_id": "test"})]
        self.assertIn("culture_interference_context_review", codes(propose(p)[0][0]))

    def test_invalid_enums_rejected(self):
        p = project()
        p["proteins"][0]["considerations"] = [dict(consideration(), status="safe")]
        with self.assertRaises(ValueError):
            validate_project(p)

    def test_prediction_cannot_be_negative_observation(self):
        p = protein()
        r = consideration()
        r.update(status="not_observed_in_context", evidence=dict(p["evidence"], kind="prediction"))
        p["considerations"] = [r]
        self.assertEqual(propose(p)[0][0]["criteria_review"][1]["status"], "unresolved")

    def test_malformed_release_context_rejected(self):
        p = project()
        p["proteins"][0]["candidate_regions"][0]["release_context"] = ["not a string"]
        with self.assertRaises(ValueError):
            validate_project(p)

    def test_approval_bound_to_context_records(self):
        p, s = protein(), scaffold()
        approval = approve(propose(p, s)[0][0])
        p["considerations"] = [consideration()]
        self.assertFalse(propose(p, s, reviews=[approval])[0][0]["fusion_sequence"])

    def test_run_exports_criteria_and_verified_manifest(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(run(project(), d)["run_dir"])
            self.assertTrue((out / "criteria_review.tsv").is_file())
            self.assertTrue(verify_bundle(out))

    def test_uniprot_processing_keeps_raw_evidence(self):
        raw = {"primaryAccession": "SYNTHETIC-001", "sequence": {"value": protein()["sequence"]}, "entryAudit": {"entryVersion": 1}, "entryType": "UniProtKB reviewed", "comments": [{"commentType": "SUBCELLULAR LOCATION", "subcellularLocations": [{"location": {"value": "Secreted"}}]}], "features": []}
        def rawfeature(kind, s, e, desc=""):
            return {"type": kind, "description": desc, "location": {"start": {"value": s}, "end": {"value": e}}}
        raw["features"] = [rawfeature("Chain", 16, 95), rawfeature("Propeptide", 11, 15)]
        p = normalize(raw)
        self.assertEqual(p["topology"], "secreted")
        self.assertEqual([f["kind"] for f in p["features"]], ["chain", "propeptide"])
        self.assertTrue(p["isoform_ambiguous"])
        raw["features"].append(rawfeature("Lipidation", 76, 76, "GPI-anchor amidated residue"))
        self.assertEqual(normalize(raw)["topology"], "gpi")


if __name__ == "__main__":
    unittest.main()
