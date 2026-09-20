"""Synthetic regression cases for the 2026-09-20 audit, not biological benchmarks."""
import copy
import tempfile
import unittest
from pathlib import Path
from test_design import approve, feature, scaffold
from test_screening import Raw, make_protein, make_set
from ecd_snipr.acquisition import disposition, select_analysis_reference
from ecd_snipr.common import read_json
from ecd_snipr.design import propose
from ecd_snipr.profile import detect_multichain_partners
from ecd_snipr.screening import _screen_entry, partition_features, run_screening
from ecd_snipr.uniprot import normalize


def raw_with_isoforms():
    raw = Raw.entry()
    raw["comments"].append({"commentType": "ALTERNATIVE PRODUCTS", "isoforms": [
        {"name": {"value": "2"}, "isoformIds": ["P00001-2"], "isoformSequenceStatus": "Displayed"},
        {"name": {"value": "1"}, "isoformIds": ["P00001-1"], "isoformSequenceStatus": "Described"}]})
    return raw


class AuditRegressions(unittest.TestCase):
    def test_displayed_isoform_is_not_assumed_to_be_one(self):
        p = select_analysis_reference(normalize(raw_with_isoforms()))
        self.assertEqual(p["isoform"], "P00001-2")
        self.assertFalse(p["isoform_ambiguous"])
        self.assertEqual(p["reference_selection"]["basis"], "database_canonical")

    def test_no_isoform_inventory_uses_explicit_canonical_marker(self):
        self.assertEqual(select_analysis_reference(normalize(Raw.entry()))["isoform"], "P00001:canonical")

    def test_noncanonical_request_cannot_relabel_canonical_coordinates(self):
        raw = raw_with_isoforms()
        raw["requested_isoform"] = "P00001-1"
        p = select_analysis_reference(normalize(raw))
        self.assertTrue(p["isoform_ambiguous"])
        self.assertEqual(p["isoform"], "P00001-1")
        self.assertIn("isoform_unresolved", {f["code"] for f in propose(p)[1]})

    def test_explicit_displayed_isoform_has_supported_coordinates(self):
        raw = raw_with_isoforms()
        raw["requested_isoform"] = "P00001-2"
        p = select_analysis_reference(normalize(raw))
        self.assertFalse(p["isoform_ambiguous"])
        self.assertEqual(p["reference_selection"]["basis"], "explicit_isoform")

    def test_foreign_isoform_features_are_not_projected(self):
        raw = raw_with_isoforms()
        foreign = copy.deepcopy(raw["features"][1])
        foreign["location"]["sequence"] = "P00001-1"
        foreign["location"].update(start={"value": 30}, end={"value": 45})
        raw["features"].append(foreign)
        p = normalize(raw)
        self.assertEqual(p["topology"], "type_i")
        self.assertEqual(len(p["features"]), 2)
        self.assertEqual(p["excluded_annotations"][0]["raw_feature"], foreign)
        raw["features"][-1]["location"]["sequence"] = "P00001-2"
        self.assertEqual(normalize(raw)["topology"], "multi_pass")

    def test_other_molecule_location_does_not_become_reference_location(self):
        raw = raw_with_isoforms()
        raw["comments"][0]["molecule"] = "Isoform 1"
        p = normalize(raw)
        self.assertEqual(p["raw_locations"], ["Cell membrane"])
        self.assertEqual(p["location"], "unknown")
        self.assertTrue(p["excluded_annotations"])
        raw["comments"][0]["molecule"] = "Isoform 2"
        self.assertEqual(normalize(raw)["location"], "plasma_membrane")

    def test_processed_fragment_subunit_is_not_applied_to_full_reference(self):
        raw = Raw.entry()
        raw["comments"].append({"commentType": "SUBUNIT", "molecule": "Released fragment",
                                "texts": [{"value": "Heterodimer with X."}]})
        p = normalize(raw)
        self.assertFalse(p.get("subunit_comments"))
        self.assertTrue(p["excluded_annotations"])

    def test_generic_membrane_is_not_organelle_exclusion(self):
        raw = Raw.entry()
        raw["features"] = []
        raw["comments"][0]["subcellularLocations"][0]["location"]["value"] = "Membrane"
        d = disposition(normalize(raw))
        self.assertEqual(d["scope_status"], "in_scope")
        self.assertEqual(d["membrane_set_membership"], "yes")
        self.assertEqual(d["natural_cell_surface_target"], "undetermined")

    def test_exact_named_main_chain_location_is_retained_but_shed_form_is_not(self):
        raw = Raw.entry()
        raw["proteinDescription"] = {"recommendedName": {"fullName": {"value": "Main protein"}}}
        raw["features"].append({"type": "Chain", "description": "Main protein",
                                "location": {"start": {"value": 11}, "end": {"value": len(raw["sequence"]["value"])}}})
        raw["comments"][0]["molecule"] = "Main protein"
        raw["comments"].append({"commentType": "SUBCELLULAR LOCATION", "molecule": "Released main protein",
                                "subcellularLocations": [{"location": {"value": "Secreted"}}]})
        p = normalize(raw)
        self.assertEqual(p["location"], "plasma_membrane")
        self.assertEqual(p["reference_locations"], ["Cell membrane"])
        self.assertEqual(len(p["excluded_annotations"]), 1)

    def test_dual_secreted_and_peripheral_location_preserves_secreted_route(self):
        raw = Raw.entry()
        raw["features"] = [{"type": "Chain", "description": "Mature chain",
                             "location": {"start": {"value": 11}, "end": {"value": 95}}}]
        raw["comments"][0]["subcellularLocations"] = [{"location": {"value": loc}}
                                                     for loc in ["Apical cell membrane", "Secreted, extracellular space"]]
        p = select_analysis_reference(normalize(raw))
        self.assertEqual((p["topology"], p["location"]), ("secreted", "plasma_membrane"))
        self.assertEqual(disposition(p)["extension_set"], "secreted_extension")
        self.assertEqual(disposition(p)["natural_cell_surface_target"], "undetermined")
        candidates, flags = propose(p)
        self.assertNotIn("secreted_location_unconfirmed", {f["code"] for f in flags})
        self.assertTrue(candidates[0]["sequence"])

    def test_named_organelle_remains_out_and_apical_membrane_is_surface(self):
        raw = Raw.entry()
        raw["features"] = []
        location = raw["comments"][0]["subcellularLocations"][0]["location"]
        location["value"] = "Mitochondrion inner membrane"
        self.assertEqual(disposition(normalize(raw))["scope_status"], "out_of_scope")
        location["value"] = "Apical cell membrane"
        self.assertEqual(normalize(raw)["location"], "plasma_membrane")

    def test_unconfirmed_analysis_reference_blocks_fusion_not_fragment(self):
        p = make_protein("SYNTHETIC-AUDIT")
        first = propose(p, scaffold())[0][0]
        approved = propose(p, scaffold(), reviews=[approve(first)])[0][0]
        self.assertTrue(approved["sequence"])
        self.assertEqual(approved["design_status"], "approved_for_assembly")
        self.assertEqual(approved["fusion_sequence"], "")
        self.assertEqual(approved["assembly_status"], "blocked")
        self.assertIn("experimental_isoform_not_confirmed", {f["code"] for f in approved["scaffold_issues"]})

    def test_confirmed_reference_still_requires_fresh_review_and_valid_scaffold(self):
        p = make_protein("SYNTHETIC-AUDIT")
        prior_review = approve(propose(p, scaffold())[0][0])
        select_analysis_reference(p, confirmation="performed")  # synthetic confirmation only
        c = propose(p, scaffold(), reviews=[prior_review])[0][0]
        self.assertFalse(c["fusion_sequence"])
        c = propose(p, scaffold(), reviews=[approve(c)])[0][0]
        self.assertTrue(c["fusion_sequence"])
        self.assertEqual(c["assembly_status"], "assembled_unvalidated")

    def test_second_uncertain_exclusion_feature_is_not_dropped(self):
        p = make_protein("SYNTHETIC-AUDIT")
        p["features"].extend([feature("propeptide", 1, 5),
                              feature("propeptide", 12, 25, boundary_status="fuzzy")])
        kept, deferred = partition_features(p)
        self.assertEqual(sum(f["kind"] == "propeptide" for f in kept), 2)
        rec, cs = _screen_entry({}, p, None, "auto")
        self.assertEqual(rec["screening_recommendation"], "insufficient_evidence")
        self.assertFalse(any(c["design_status"] != "blocked" for c in cs))

    def test_blocked_alternative_cannot_inherit_positive_primary(self):
        p = make_protein("SYNTHETIC-AUDIT")
        p["candidate_regions"] = [dict(feature("region", 11, 82), rationale="Synthetic TM-overlapping alternative")]
        with tempfile.TemporaryDirectory() as tmp:
            make_set(Path(tmp) / "set", [p])
            result = run_screening(Path(tmp) / "set", Path(tmp) / "out", pilot=True)
            candidates = read_json(Path(result["run_dir"]) / "candidates.json")
        bad = next(c for c in candidates if c["design_status"] == "blocked")
        self.assertEqual(bad["screening_recommendation"], "no_standard_route")
        self.assertEqual(bad["protein_screening_recommendation"], "standard_candidate")
        self.assertIn("retained_transmembrane", bad["screening_reason_codes"])

    def test_domain_alternative_is_evaluated_separately_from_full_ecd(self):
        p = make_protein("SYNTHETIC-AUDIT")
        rec, cs = _screen_entry({}, p, None, "full_ecd_and_domains")
        self.assertEqual(rec["screening_recommendation"], "standard_candidate")
        self.assertTrue(all(c["screening_recommendation"] == "conditional_candidate"
                            for c in cs if c["origin"] == "domain_alternative"))

    def test_negative_heteromer_sentence_is_not_positive_evidence(self):
        self.assertEqual(detect_multichain_partners([{"text": "Does not form heterodimers."}]), [])
        self.assertEqual(detect_multichain_partners([{"text": "Forms a heterodimer with X."}]),
                         ["Forms a heterodimer with X."])

    def test_heteromer_is_context_not_proven_ecd_dependency(self):
        p = make_protein("SYNTHETIC-AUDIT")
        p["subunit_comments"] = [{"text": "Heterodimer with X."}]
        c = propose(p)[0][0]
        codes = {r["code"] for r in c["risks"]}
        self.assertIn("native_heteromer_context_requires_review", codes)
        self.assertNotIn("multichain_partner_required", codes)


if __name__ == "__main__":
    unittest.main()
