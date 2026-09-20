import copy
import tempfile
import unittest
from pathlib import Path
from test_design import feature
from test_screening import make_protein, make_set
from ecd_snipr.common import read_json
from ecd_snipr.harness import verify_bundle
from ecd_snipr.methodology import evidence_support, candidate_tradeoff, receiver_review, review_queue
from ecd_snipr.screening import pick_primary, run_screening


class EvidenceTests(unittest.TestCase):
    def test_reviewed_unspecified_not_experimental(self):
        self.assertEqual(evidence_support({"kind": "curated_annotation"})["support_types"], ["curated_support_unspecified"])

    def test_eco_labels_not_entry_labels(self):
        for code, label in [("ECO:0000269", "experimental_annotation"),
                            ("ECO:0000250", "similarity_transfer"), ("ECO:0000255", "sequence_model_inference")]:
            evidence = {"kind": "curated_annotation", "eco": [{"evidenceCode": code, "source": "123"}]}
            r = evidence_support(evidence)
            self.assertEqual(r["support_types"], [label])
            self.assertEqual(r["provenance"], evidence)

    def test_unknown_code_preserved_not_promoted(self):
        r = evidence_support({"kind": "curated_annotation", "eco": ["ECO:9999999"]})
        self.assertEqual(r["support_types"], ["unclassified_eco"])
        self.assertEqual(r["eco_codes"], ["ECO:9999999"])

    def test_mixed_evidence_not_flattened(self):
        r = evidence_support({"eco": ["ECO:0000269", "ECO:0000255"]})
        self.assertEqual(len(r["support_types"]), 2)

    def test_fixture_not_promoted_by_eco(self):
        r = evidence_support({"kind": "synthetic_fixture", "eco": ["ECO:0000269"]})
        self.assertEqual(r["support_types"], ["synthetic_only"])

    def test_optional_eco_shapes(self):
        for value in ("ECO:0000250", {"evidenceCode": "ECO:0000250"}):
            self.assertEqual(evidence_support({"eco": value})["support_types"], ["similarity_transfer"])
        self.assertEqual(evidence_support({"eco": None})["support_types"], ["missing"])


def cand(cid="a", form="full_ecd", codes=()):
    return {"candidate_id": cid, "sequence": "AAA", "design_status": "needs_review", "antigen_form_type": form,
            "risks": [{"code": code, "severity": "review"} for code in codes]}


class ComparisonTests(unittest.TestCase):
    def test_full_form_preferred_not_warning_counts(self):
        full = cand(codes=["long_fragment_geometry_review", "potential_n_glycosylation"])
        alt = cand("b", "domain_fragment")
        p, a = pick_primary([alt, full])
        self.assertEqual(p["candidate_id"], "a")
        self.assertEqual(a[0]["deciding_axis"], "antigen_form_priority")

    def test_integrity_disruption_changes_selection(self):
        full = cand(codes=["domain_cut"])
        alt = cand("b", "domain_fragment")
        p, a = pick_primary([full, alt])
        self.assertEqual(p["candidate_id"], "b")
        self.assertEqual(a[0]["deciding_axis"], "annotated_integrity_disruption")

    def test_known_loss_not_ignored(self):
        a, b = cand(), cand("b")
        a["epitope_review"] = [{"name": "mapped", "state": "removed"}]
        p, _ = pick_primary([a, b])
        self.assertEqual(p["candidate_id"], "b")
        self.assertEqual(candidate_tradeoff(b)["mapped_epitope_status"], "unknown")

    def test_blocked_cannot_win(self):
        a, b = cand(), cand("b", "domain_fragment")
        a["design_status"] = "blocked"
        self.assertEqual(pick_primary([a, b])[0]["candidate_id"], "b")

    def test_no_warning_count_tiebreak(self):
        a, b = cand(codes=["extra_warning"]), cand("b")
        p, alts = pick_primary([b, a])
        self.assertEqual(p["candidate_id"], "a")
        self.assertEqual(alts[0]["deciding_axis"], "candidate_id_tiebreak")

    def test_receiver_four_unmeasured_endpoints(self):
        r = receiver_review(cand())
        self.assertEqual(len(r["endpoints"]), 4)
        self.assertEqual(set(r["endpoints"].values()), {"not_measured"})
        self.assertFalse(r["fusion_sequence_authorized"])
        self.assertEqual(r["checks"][0]["status"], "not_evaluated")


class ScreeningMethodologyTests(unittest.TestCase):
    def screen(self, proteins, **kwargs):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        make_set(root / "set", proteins)
        r = run_screening(root / "set", root / "out", pilot=True, **kwargs)
        run = Path(r["run_dir"])
        return r, read_json(run / "screening.json"), read_json(run / "candidates.json")

    def test_no_scaffold_no_review_end_to_end(self):
        result, rs, cs = self.screen([make_protein("P98001")])
        self.assertEqual(rs[0]["screening_recommendation"], "standard_candidate")
        self.assertEqual(rs[0]["route_diagnostic"]["state"], "supported_candidate_unvalidated")
        self.assertTrue(verify_bundle(result["run_dir"]))
        for filename in ("candidate_comparison.tsv", "receiver_review_plan.tsv", "manual_review_queue.tsv"):
            self.assertTrue((Path(result["run_dir"]) / filename).stat().st_size)
        self.assertFalse(any(c["fusion_sequence"] for c in cs))

    def test_feature_support_not_entry_review_status(self):
        p = make_protein("P98002")
        p["evidence"]["kind"] = "curated_annotation"
        p["features"][1]["evidence"] = {"kind": "prediction", "source": "synthetic://test", "version": "1", "eco": ["ECO:0000255"]}
        _, rs, cs = self.screen([p])
        self.assertEqual(rs[0]["screening_recommendation"], "standard_candidate")
        self.assertEqual(cs[0]["annotation_support"]["boundary"]["support_types"], ["sequence_model_inference"])

    def test_missing_external_annotation_is_insufficient(self):
        p = make_protein("P98003")
        p["features"] = [f for f in p["features"] if f["kind"] != "extracellular"]
        _, rs, _ = self.screen([p])
        self.assertEqual(rs[0]["screening_recommendation"], "insufficient_evidence")
        self.assertIn("extracellular_boundary_missing", rs[0]["reason_codes"])

    def test_multipass_missing_domain_is_route_gap_not_impossibility(self):
        p = make_protein("P98004")
        p["topology"] = "multi_pass"
        p["features"] = [feature("extracellular", 10, 20), feature("transmembrane", 21, 39),
                         feature("extracellular", 40, 50), feature("transmembrane", 51, 60)]
        _, rs, cs = self.screen([p])
        self.assertFalse(cs)
        self.assertEqual(rs[0]["route_diagnostic"]["state"], "no_supported_route_in_current_annotations")
        self.assertFalse(rs[0]["route_diagnostic"]["biological_impossibility_claim"])

    def test_domain_cut_and_disulfide_are_distinct_conditional_reasons(self):
        p = make_protein("P98005")
        p["features"].append(feature("disulfide", 30, 80))
        _, rs, cs = self.screen([p])
        self.assertIn("annotated_disulfide_crosses_boundary", rs[0]["reason_codes"])
        self.assertEqual(rs[0]["screening_recommendation"], "conditional_candidate")
        self.assertIn("disulfide_partner_removed", cs[0]["candidate_comparison"]["integrity_disruptions"])

    def test_unknown_epitope_does_not_remove_candidate(self):
        p = make_protein("P98006")
        p["features"] = [f for f in p["features"] if f["kind"] != "epitope"]
        result, rs, cs = self.screen([p])
        self.assertEqual(rs[0]["screening_recommendation"], "standard_candidate")
        self.assertEqual(cs[0]["candidate_comparison"]["mapped_epitope_status"], "unknown")
        self.assertEqual(result["summary"]["methodology"]["primary_with_unknown_epitopes"], 1)

    def test_invalid_sequence_not_rescued_by_methodology(self):
        p = make_protein("P98007")
        p["sequence"] = "X*"
        _, rs, cs = self.screen([p])
        self.assertFalse(cs)
        self.assertEqual(rs[0]["screening_recommendation"], "insufficient_evidence")

    def test_queue_limits_determinism_and_empty_labels(self):
        _, rs, cs = self.screen([make_protein("P98008"), make_protein("P98009")])
        queue = review_queue(rs, cs, 1)
        self.assertEqual(queue, review_queue(list(reversed(rs)), list(reversed(cs)), 1))
        self.assertEqual(len(queue), 1)
        self.assertEqual(queue[0]["review_status"], "pending")
        self.assertEqual(queue[0]["induced_response"], "")
        self.assertEqual(queue[0]["stratum_population"], 2)

    def test_queue_range_validation(self):
        for value in (0, -1, 101, True):
            with self.assertRaises(ValueError):
                review_queue([], [], value)

    def test_review_options_part_of_run_identity(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        make_set(root / "set", [make_protein("P98010")])
        a = run_screening(root / "set", root / "out", review_per_stratum=1)
        b = run_screening(root / "set", root / "out", review_per_stratum=2)
        self.assertNotEqual(a["run_dir"], b["run_dir"])


if __name__ == "__main__":
    unittest.main()
