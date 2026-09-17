import copy
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from ecd_snipr.common import read_json, translate
from ecd_snipr.design import audit_construct, propose, scaffold_check


def project():
    return read_json(ROOT / "examples/synthetic_project.json")


def protein():
    return project()["proteins"][0]


def feature(kind, start, end, **extra):
    return dict(kind=kind, start=start, end=end, evidence=protein()["evidence"], **extra)


def scaffold():
    roles = [("signal_peptide", "MAAAA"), ("antigen_slot", ""), ("linker", "GSG"), ("transmembrane", "LLLLL"), ("juxtamembrane", "AAA"), ("transcription_factor", "RRR")]
    return {"scaffold_id": "SYNTHETIC-NONFUNCTIONAL", "version": "1", "source": "synthetic://test", "host_cell": "synthetic-host", "reporter": "synthetic-reporter", "status": "synthetic_fixture", "signal_policy": "scaffold_signal_only", "junctions_verified": True,
            "configuration": "antibody_sender_antigen_receiver",
            "modules": [{"role": role, "name": role, "sequence": seq, "source": "synthetic://test"} for role, seq in roles]}


def approve(candidate):
    return {"review_key": candidate["review_key"], "decision": "approve", "reviewer": "synthetic-test-only", "rationale": "Synthetic test approval, not experimental approval", "reviewed_at": "2026-01-01T00:00:00Z", "acknowledged_risks": [r["code"] for r in candidate["risks"] if r["severity"] == "review"]}


def codes(candidate):
    return {r["code"] for r in candidate["risks"]}


class DesignTests(unittest.TestCase):
    def test_baseline_and_complete_domain_alternatives(self):
        p = protein()
        candidates, _ = propose(p)
        self.assertEqual([(c["start"], c["end"]) for c in candidates], [(11, 76), (15, 40), (45, 70)])
        self.assertEqual(candidates[0]["sequence"], p["sequence"][10:76])
        self.assertEqual(candidates[0]["domains_retained"], ["synthetic-domain-A", "synthetic-domain-B"])

    def test_no_scaffold_no_fusion(self):
        candidates, _ = propose(protein())
        self.assertTrue(all(not c["fusion_sequence"] for c in candidates))
        self.assertEqual(candidates[0]["scaffold_issues"][0]["code"], "scaffold_missing")

    def test_known_epitope_loss(self):
        candidates, _ = propose(protein())
        self.assertIn("known_epitope_loss", codes(candidates[1]))
        self.assertNotIn("known_epitope_loss", codes(candidates[0]))

    def test_unknown_epitope_not_absent(self):
        p = protein()
        p["features"] = [f for f in p["features"] if f["kind"] != "epitope"]
        self.assertIn("epitope_coverage_unknown", codes(propose(p)[0][0]))

    def test_domain_cut_and_disulfide(self):
        p = protein()
        p["features"].append(feature("disulfide", 22, 55))
        p["candidate_regions"] = [dict(feature("region", 20, 30), rationale="Sourced synthetic truncation")]
        candidate = propose(p)[0][-1]
        self.assertIn("domain_cut", codes(candidate))
        self.assertIn("disulfide_partner_removed", codes(candidate))

    def test_tm_overlap_blocks(self):
        p = protein()
        p["candidate_regions"] = [dict(feature("region", 11, 82), rationale="Deliberately invalid test")]
        self.assertEqual(propose(p)[0][-1]["design_status"], "blocked")
        self.assertIn("retained_transmembrane", codes(propose(p)[0][-1]))

    def test_signal_overlap_blocks(self):
        p = protein()
        p["features"][1]["start"] = 5
        self.assertIn("retained_signal_peptide", codes(propose(p)[0][0]))

    def test_fuzzy_not_rounded(self):
        p = protein()
        p["candidate_regions"][0]["boundary_status"] = "fuzzy"
        c = propose(p)[0][-1]
        self.assertEqual(c["sequence"], "")
        self.assertEqual(c["design_status"], "blocked")

    def test_out_of_range(self):
        p = protein()
        p["candidate_regions"] = [dict(feature("region", 2, 500), rationale="Test")]
        self.assertEqual(propose(p)[0][-1]["sequence"], "")

    def test_unknown_isoform_blocks_assembly(self):
        p = protein()
        p["isoform_ambiguous"] = True
        c = propose(p, scaffold())[0][0]
        self.assertIn("isoform_unresolved", codes(c))
        self.assertEqual(propose(p, scaffold(), reviews=[approve(c)])[0][0]["fusion_sequence"], "")

    def test_multipass_no_loop_concatenation(self):
        p = protein()
        p["topology"] = "multi_pass"
        p["features"] = [feature("extracellular", 10, 20), feature("extracellular", 40, 50), feature("transmembrane", 21, 39)]
        p["candidate_regions"] = []
        candidates, issues = propose(p)
        self.assertEqual(candidates, [])
        self.assertIn("no_supported_continuous_candidate", {f["code"] for f in issues})

    def test_multipass_domain_requires_review(self):
        p = protein()
        p["topology"] = "multi_pass"
        self.assertTrue(all(c["origin"] != "baseline" for c in propose(p)[0]))
        self.assertIn("discontinuous_or_multipass", codes(propose(p)[0][0]))

    def test_type_ii_orientation(self):
        p = protein()
        p["topology"] = "type_ii"
        self.assertIn("type_ii_attachment_orientation_change", codes(propose(p)[0][0]))

    def test_contradictory_topology_blocks(self):
        p = protein()
        p["topology"] = "type_ii"
        self.assertIn("topology_orientation_conflict", codes(propose(p)[0][0]))
        self.assertEqual(propose(p)[0][0]["design_status"], "blocked")

    def test_true_type_ii_extracts_C_terminal_ecd(self):
        p = protein()
        p["topology"] = "type_ii"
        p["features"] = [feature("cytoplasmic", 1, 10), feature("transmembrane", 11, 20), feature("extracellular", 21, 95), feature("domain", 30, 80)]
        p["candidate_regions"] = []
        candidates, _ = propose(p)
        self.assertEqual(candidates[0]["sequence"], p["sequence"][20:95])
        self.assertNotEqual(candidates[0]["design_status"], "blocked")

    def test_gpi_signal_overlap_blocks(self):
        p = protein()
        p["topology"] = "gpi"
        p["features"].append(feature("gpi_signal", 70, 95))
        self.assertIn("retained_gpi_signal", codes(propose(p)[0][0]))

    def test_gpi_review(self):
        p = protein()
        p["topology"] = "gpi"
        self.assertIn("gpi_anchor_replaced", codes(propose(p)[0][0]))

    def test_lumen_not_surface(self):
        p = protein()
        p["location"] = "other_membrane"
        p["features"][1]["kind"] = "lumenal"
        p["candidate_regions"] = []
        candidates, _ = propose(p)
        self.assertEqual(candidates, [])

    def test_short_long_cysteine_and_glycosylation_rules(self):
        p = protein()
        p["sequence"] = "NAT" + "C" * 92
        candidates, _ = propose(p, rules={"long_length": 60})
        self.assertIn("long_fragment_geometry_review", codes(candidates[0]))
        self.assertIn("cysteine_rich_folding_review", codes(candidates[0]))
        p["sequence"] = p["sequence"][:20] + "NAT" + p["sequence"][23:]
        self.assertIn("potential_n_glycosylation", codes(propose(p)[0][0]))

    def test_prediction_not_experimental(self):
        p = protein()
        p["features"][1]["evidence"]["kind"] = "prediction"
        self.assertIn("prediction_requires_annotation_review", codes(propose(p)[0][0]))

    def test_missing_provenance_blocks_without_crash(self):
        p = protein()
        del p["evidence"]
        self.assertEqual(propose(p)[0][0]["design_status"], "blocked")

    def test_source_free_region_blocks(self):
        p = protein()
        p["candidate_regions"] = [{"start": 12, "end": 66, "rationale": "LLM guess"}]
        self.assertIn("candidate_support_missing", codes(propose(p)[0][-1]))

    def test_scaffold_validation(self):
        s = scaffold()
        self.assertEqual(scaffold_check(s), [])
        s["modules"].insert(0, s["modules"][0])
        self.assertTrue(scaffold_check(s))

    def test_wrong_configuration(self):
        s = scaffold()
        s["configuration"] = "antigen_sender"
        self.assertTrue(scaffold_check(s))

    def test_public_reference_is_not_actual_scaffold(self):
        reference = read_json(ROOT / "references/public_reference_architecture.json")
        candidates, _ = propose(protein(), reference)
        self.assertEqual(len(candidates), 3)
        self.assertNotEqual(candidates[0]["design_status"], "blocked")
        self.assertFalse(candidates[0]["fusion_sequence"])
        self.assertIn("reference_architecture_not_lab_scaffold", {f["code"] for f in candidates[0]["scaffold_issues"]})

    def test_public_reference_role_and_sequence_exclusion(self):
        reference = read_json(ROOT / "references/public_reference_architecture.json")
        self.assertFalse(reference["assembly_authorized"])
        self.assertEqual(reference["current_lab_equivalence"], "unconfirmed")
        sender = next(c for c in reference["components"] if c["component"] == "excluded_sender_sequences")
        self.assertFalse(sender["allowed_as_receiver_scaffold"])
        self.assertTrue(all(s["url"].startswith("https://") for s in reference["sources"]))

    def test_review_required_and_correct_assembly(self):
        p, s = protein(), scaffold()
        c = propose(p, s)[0][0]
        self.assertEqual(c["assembly_status"], "awaiting_review")
        c = propose(p, s, reviews=[approve(c)])[0][0]
        self.assertEqual(c["fusion_sequence"], "MAAAA" + p["sequence"][10:76] + "GSGLLLLLAAARRR")
        self.assertTrue(c["synthetic_only"])

    def test_review_invalidated_by_scaffold_change(self):
        p, s = protein(), scaffold()
        review = approve(propose(p, s)[0][0])
        s["version"] = "2"
        self.assertFalse(propose(p, s, reviews=[review])[0][0]["fusion_sequence"])

    def test_unacknowledged_risk_blocks_approval(self):
        p, s = protein(), scaffold()
        c = propose(p, s)[0][1]
        review = approve(c)
        review["acknowledged_risks"] = []
        self.assertFalse(propose(p, s, reviews=[review])[0][1]["fusion_sequence"])

    def test_translation_and_internal_stop(self):
        self.assertEqual(translate("ATGGCTTAA"), "MA")
        with self.assertRaises(ValueError):
            translate("ATGTAGGCT")

    def test_historical_sequence_mismatch(self):
        p = protein()
        record = {"protein_id": p["protein_id"], "construct_id": "TEST", "start": 11, "end": 76, "antigen_sequence": "A" * 66}
        result = audit_construct(record, {p["protein_id"]: p})
        self.assertIn("actual_sequence_differs_from_reference_interval", result["audit_issues"])

    def test_unknown_residue_not_fabricated(self):
        p = protein()
        p["sequence"] += "X"
        self.assertEqual(propose(p)[0], [])


if __name__ == "__main__":
    unittest.main()
