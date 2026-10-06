"""Synthetic counterexamples for biological-unit scope; no functional labels."""
import copy
import re
import tempfile
import unittest
from pathlib import Path
from test_design import feature, approve, scaffold
from test_screening import Raw, make_protein, make_set
from ecd_snipr.acquisition import select_analysis_reference
from ecd_snipr.antigen_context import assess_antigen_context
from ecd_snipr.common import read_json
from ecd_snipr.design import propose
from ecd_snipr.harness import verify_bundle
from ecd_snipr.screening import _screen_entry, partition_features, recommend, run_screening
from ecd_snipr.uniprot import normalize


def chain_raw(name, start, end, sequence=None):
    value = {"type": "Chain", "description": name,
             "location": {"start": {"value": start}, "end": {"value": end}}}
    if sequence:
        value["location"]["sequence"] = sequence
    return value


def unit_raw():
    raw = Raw.entry()
    raw["features"].append(chain_raw("Secreted alpha", 11, 30))
    raw["comments"].append({"commentType": "SUBCELLULAR LOCATION", "molecule": "Secreted alpha",
                            "subcellularLocations": [{"location": {"value": "Secreted, extracellular space"}}]})
    return raw


class ProductMappingTests(unittest.TestCase):
    def test_product_localization_not_promoted_to_precursor(self):
        p = normalize(unit_raw())
        self.assertEqual(p["location"], "plasma_membrane")
        self.assertFalse(p["has_secreted_reference_annotation"])
        self.assertEqual(p["processed_products"][0]["scope"], "processed_product_not_whole_reference")
        self.assertEqual(p["processed_products"][0]["start"], 11)

    def test_foreign_isoform_chain_does_not_acquire_context(self):
        raw = unit_raw()
        raw["features"][-1]["location"]["sequence"] = "P00001-9"
        self.assertFalse(normalize(raw)["processed_products"])

    def test_fuzzy_chain_cannot_define_product(self):
        raw = unit_raw()
        raw["features"][-1]["location"]["end"]["modifier"] = "UNCERTAIN"
        self.assertFalse(normalize(raw)["processed_products"])

    def test_ambiguous_product_name_not_guessed(self):
        raw = unit_raw()
        raw["features"].append(chain_raw("Secreted alpha", 11, 45))
        self.assertFalse(normalize(raw)["processed_products"])

    def test_no_name_match_not_guessed(self):
        raw = unit_raw()
        raw["comments"][-1]["molecule"] = "Another chain"
        self.assertFalse(normalize(raw)["processed_products"])

    def test_exact_and_fuzzy_duplicate_name_still_ambiguous(self):
        raw = unit_raw()
        duplicate = chain_raw("Secreted alpha", 11, 45)
        duplicate["location"]["end"]["modifier"] = "UNCERTAIN"
        raw["features"].append(duplicate)
        self.assertFalse(normalize(raw)["processed_products"])

    def test_one_tm_multiple_external_segments_orientation_known(self):
        raw = Raw.entry(ext=(11, 30))
        raw["features"].append({"type": "Topological domain", "description": "Extracellular",
                                "location": {"start": {"value": 35}, "end": {"value": 70}}})
        p = select_analysis_reference(normalize(raw))
        self.assertEqual(p["topology"], "type_i")
        # Known orientation is not permission to stitch segments.
        cs, _ = propose(p)
        self.assertFalse(cs)

    def test_opposite_sides_not_assigned_single_orientation(self):
        raw = Raw.entry()
        raw["features"].append({"type": "Topological domain", "description": "Extracellular",
                                "location": {"start": {"value": 91}, "end": {"value": 95}}})
        self.assertEqual(normalize(raw)["topology"], "unknown")

    def test_repeat_import_is_not_domain_route(self):
        raw = Raw.entry()
        raw["features"].append({"type": "Repeat", "description": "Synthetic repeat",
                                "location": {"start": {"value": 12}, "end": {"value": 32}}})
        p = select_analysis_reference(normalize(raw))
        p["candidate_policy"] = "full_ecd_and_domains"
        self.assertEqual(p["features"][-1]["kind"], "repeat")
        self.assertFalse(any(c["antigen_form_type"] == "domain_fragment" for c in propose(p)[0]))


class ContextTests(unittest.TestCase):
    def test_nested_chains_not_processing_split(self):
        p = make_protein("P98050")
        p["features"] += [feature("chain", 1, 95), feature("chain", 11, 95)]
        context, flags = assess_antigen_context(p, (11, 70))
        self.assertFalse(context["processed_chain_pairs_spanned"])
        self.assertNotIn("processed_chain_segments_spanned", {f["code"] for f in flags})

    def test_disjoint_chain_segments_require_processing_review(self):
        p = make_protein("P98051")
        p["features"] += [feature("chain", 11, 45), feature("chain", 46, 95)]
        cs, fs = propose(p)
        rec = recommend(p, cs, fs)
        self.assertEqual(rec["screening_recommendation"], "conditional_candidate")
        self.assertIn("processed_chain_segments_spanned", rec["reason_codes"])
        self.assertTrue(cs[0]["antigen_context"]["processed_chain_pairs_spanned"])

    def test_external_product_omitted_is_not_all_gene_epitopes(self):
        p = select_analysis_reference(normalize(unit_raw()))
        context, flags = assess_antigen_context(p, (35, 70))
        self.assertEqual(context["extracellular_mature_products"][0]["coverage"], "omitted")
        self.assertIn("extracellular_processed_product_not_covered", {f["code"] for f in flags})

    def test_product_sequence_contained_not_binding_proof(self):
        p = select_analysis_reference(normalize(unit_raw()))
        context, flags = assess_antigen_context(p, (11, 70))
        self.assertEqual(context["extracellular_mature_products"][0]["coverage"], "sequence_contained")
        self.assertEqual(context["recognition_retention"], "not_measured")
        self.assertNotIn("extracellular_processed_product_not_covered", {f["code"] for f in flags})

    def test_product_crossing_tm_not_claimed_soluble_unit(self):
        raw = unit_raw()
        raw["features"][-1]["location"]["end"]["value"] = 95
        context, _ = assess_antigen_context(normalize(raw), (11, 70))
        self.assertFalse(context["extracellular_mature_products"])

    def test_repeat_cut_not_failure_or_autonomous_domain(self):
        p = make_protein("P98052")
        p["features"].append(feature("repeat", 60, 80))
        cs, fs = propose(p)
        self.assertNotEqual(cs[0]["design_status"], "blocked")
        self.assertEqual(recommend(p, cs, fs)["screening_recommendation"], "conditional_candidate")
        self.assertEqual(cs[0]["antigen_context"]["overlapping_repeats"][0]["coverage"], "cut")

    def test_repeat_selection_reason_names_actual_tradeoff(self):
        p = make_protein("P98057")
        p["features"] = [f for f in p["features"] if f["kind"] not in {"domain", "epitope", "disulfide"}]
        p["features"] += [feature("repeat", 60, 80), feature("domain", 12, 29)]
        p["candidate_policy"] = "full_ecd_and_domains"
        cs, fs = propose(p)
        rec = recommend(p, cs, fs)
        chosen = next(c for c in cs if c["candidate_id"] == rec["primary_candidate_id"])
        self.assertEqual((chosen["start"], chosen["end"]), (12, 29))
        self.assertIn("重复单元", rec["alternates"][0]["reason_not_primary"])
        self.assertEqual(chosen["antigen_context"]["recognition_retention"], "not_measured")

    def test_missing_epitope_source_is_not_negative_search(self):
        p = make_protein("P98053")
        p["features"] = [f for f in p["features"] if f["kind"] != "epitope"]
        context, _ = assess_antigen_context(p, (11, 70))
        self.assertEqual(context["epitope_evidence_scope"], "no_mapped_annotations_supplied")
        self.assertEqual(context["external_epitope_database_search"], "not_performed")

    def test_reason_codes_do_not_alias_missing_info(self):
        p = make_protein("P98054")
        rec = recommend(p, [], [{"code": "topology_unknown", "severity": "block"}])
        rec["missing_info"].append("A prose note")
        self.assertNotIn("A prose note", rec["reason_codes"])

    def test_review_flags_still_gate_assembly(self):
        p = make_protein("P98055")
        p["features"] += [feature("chain", 11, 45), feature("chain", 46, 95)]
        p["reference_selection"]["experimental_isoform_confirmation"] = "performed"
        candidates, _ = propose(p, scaffold())
        self.assertFalse(any(c["fusion_sequence"] for c in candidates))

    def test_output_context_is_manifest_verified(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            p = make_protein("P98056")
            p["features"] += [feature("chain", 11, 45), feature("chain", 46, 95)]
            make_set(root / "set", [p])
            result = run_screening(root / "set", root / "out", pilot=True)
            out = Path(result["run_dir"])
            self.assertTrue(verify_bundle(out))
            for name in ("antigen_context.tsv", "processed_product_review.tsv"):
                self.assertIn(name, read_json(out / "manifest.json")["artifacts"])
            self.assertEqual(result["summary"]["methodology"]["primary_with_processing_or_repertoire_flags"], 1)
            for r in read_json(out / "screening.json"):
                self.assertTrue(all(re.fullmatch(r"[a-z0-9_]+", c) for c in r["reason_codes"]))


class IntramembraneTests(unittest.TestCase):
    def test_embedded_segment_is_not_second_transmembrane_span(self):
        raw = Raw.entry(tm=(21, 40), ext=(41, 65))
        raw["features"] += [
            {"type": "Intramembrane", "location": {"start": {"value": 66}, "end": {"value": 80}}},
            {"type": "Topological domain", "description": "Extracellular",
             "location": {"start": {"value": 81}, "end": {"value": 95}}}]
        p = select_analysis_reference(normalize(raw))
        self.assertEqual(p["topology"], "type_ii")
        self.assertEqual(sum(f["kind"] == "transmembrane" for f in p["features"]), 1)
        cs, fs = propose(p)
        self.assertFalse(cs)  # No joining of extracellular segments across embedded sequence.
        r = recommend(p, cs, fs)
        self.assertEqual(r["screening_recommendation"], "no_standard_route")
        self.assertIn("intramembrane_segmented_route_requires_review", r["reason_codes"])
        self.assertNotIn("segmented_single_pass_needs_product_specific_design", r["reason_codes"])

    def test_embedded_sequence_overlap_blocks_soluble_candidate(self):
        p = make_protein("P98101")
        p["features"].append(feature("intramembrane", 45, 55))
        cs, _ = propose(p)
        self.assertTrue(cs)
        self.assertTrue(all(c["design_status"] == "blocked" for c in cs))
        self.assertIn("retained_intramembrane", {r["code"] for r in cs[0]["risks"]})

    def test_fuzzy_embedded_exclusion_is_not_deferred(self):
        p = make_protein("P98102")
        f = feature("intramembrane", 45, 55)
        f["boundary_status"] = "fuzzy"
        p["features"].append(f)
        engine, deferred = partition_features(p)
        self.assertIn(f, engine)
        self.assertFalse(deferred)
        cs, fs = propose(p)
        self.assertEqual(recommend(p, cs, fs)["screening_recommendation"], "insufficient_evidence")

    def test_external_domain_does_not_require_removing_safe_route(self):
        p = make_protein("P98103")
        p["features"] = [f for f in p["features"] if f["kind"] not in {"domain", "extracellular"}]
        p["features"] += [feature("extracellular", 11, 30), feature("intramembrane", 31, 45),
                          feature("extracellular", 46, 70), feature("domain", 12, 29)]
        p["candidate_policy"] = "full_ecd_and_domains"
        cs, fs = propose(p)
        self.assertEqual([(c["start"], c["end"]) for c in cs], [(12, 29)])
        self.assertNotEqual(cs[0]["design_status"], "blocked")
        self.assertEqual(recommend(p, cs, fs)["screening_recommendation"], "conditional_candidate")

    def test_embedded_product_not_claimed_separable_soluble_unit(self):
        raw = unit_raw()
        raw["features"].append({"type": "Intramembrane", "location": {
            "start": {"value": 20}, "end": {"value": 25}}})
        context, _ = assess_antigen_context(normalize(raw), (35, 70))
        self.assertFalse(context["extracellular_mature_products"])
        self.assertEqual(context["reference_intramembrane_segments"][0]["start"], 20)
