"""v0.4.0 deep-screening tests.

Synthetic entries exercising each new deterministic check: fuller UniProt
feature extraction, per-candidate molecular profile (cysteine parity, disulfide
containment, glycosylation, domain cuts, site/variant overlap), boundary
precision analysis (exact / tolerance / conflict / missing), multichain
SUBUNIT evidence, multi-pass loop enumeration, and the scoped
prediction-evidence review. Warning-level codes must never change a class.
"""
import copy
import tempfile
import unittest
from pathlib import Path
from test_design import feature, protein
from test_screening import Raw, make_protein, make_set
from ecd_snipr.common import read_json
from ecd_snipr.design import propose
from ecd_snipr.profile import (boundary_analysis, detect_multichain_partners,
                               enumerate_external_loops, molecular_profile, n_glyco_sequons)
from ecd_snipr.screening import run_screening
from ecd_snipr.uniprot import normalize


def codes(candidate):
    return {r["code"] for r in candidate["risks"]}


def screen_one(p, **kwargs):
    temp = tempfile.TemporaryDirectory()
    setdir = Path(temp.name) / "set"
    make_set(setdir, [p])
    result = run_screening(setdir, Path(temp.name) / "screen", **kwargs)
    records = read_json(Path(result["run_dir"]) / "screening.json")
    candidates = read_json(Path(result["run_dir"]) / "candidates.json")
    return temp, result, records, candidates  # caller must keep temp alive


class FullerExtractionTests(unittest.TestCase):
    def test_normalize_parses_all_feature_kinds_with_evidence(self):
        raw = Raw.entry("P90001", "DEEP")
        raw["annotationScore"] = 4.0
        raw["features"] += [
            {"type": "Region", "description": "Interaction with X",
             "location": {"start": {"value": 15}, "end": {"value": 25}}},
            {"type": "Motif", "description": "NLS-like",
             "location": {"start": {"value": 30}, "end": {"value": 34}}},
            {"type": "Glycosylation", "description": "N-linked (GlcNAc...) asparagine",
             "location": {"start": {"value": 40}, "end": {"value": 40}}, "featureId": "CARB_1",
             "evidences": [{"evidenceCode": "ECO:0000256"}]},
            {"type": "Site", "description": "Cleavage; by thrombin",
             "location": {"start": {"value": 50}, "end": {"value": 51}}},
            {"type": "Binding site", "description": "Substrate",
             "location": {"start": {"value": 55}, "end": {"value": 58}}},
            {"type": "Active site", "description": "Proton acceptor",
             "location": {"start": {"value": 60}, "end": {"value": 60}}},
            {"type": "Natural variant", "description": "R -> H (in a population)",
             "location": {"start": {"value": 45}, "end": {"value": 45}}, "featureId": "VAR_1"},
            {"type": "Mutagenesis", "description": "A->G: loss of binding",
             "location": {"start": {"value": 48}, "end": {"value": 48}}, "featureId": "MUT_1"},
            {"type": "Lipidation", "description": "S-palmitoyl cysteine",
             "location": {"start": {"value": 65}, "end": {"value": 65}}},
        ]
        p = normalize(raw)
        kinds = {f["kind"] for f in p["features"]}
        self.assertTrue({"region", "motif", "glycosylation_site", "site", "binding_site",
                         "active_site", "variant", "mutagenesis", "lipidation"} <= kinds)
        gly = next(f for f in p["features"] if f["kind"] == "glycosylation_site")
        self.assertEqual((gly["start"], gly["end"]), (40, 40))
        self.assertEqual(gly["name"], "N-linked (GlcNAc...) asparagine")
        self.assertEqual(gly["feature_id"], "CARB_1")
        self.assertEqual(gly["evidence"]["kind"], "prediction")  # ECO:0000256
        var = next(f for f in p["features"] if f["kind"] == "variant")
        self.assertEqual(var["evidence"]["kind"], "curated_annotation")  # reviewed entry, no prediction ECO
        self.assertEqual(p["annotation_score"], 4.0)
        self.assertTrue(p["reviewed"])
        # A non-GPI lipidation is informational, not a gpi attachment site.
        self.assertFalse(any(f["kind"] == "gpi_attachment_site" for f in p["features"]))

    def test_normalize_records_subunit_function_ptm_comments(self):
        raw = Raw.entry("P90002", "DEEP2")
        raw["comments"] += [
            {"commentType": "SUBUNIT", "texts": [{"value": "Monomer. Forms a heterodimer with PARTNERX."}],
             "evidences": [{"evidenceCode": "ECO:0000269"}]},
            {"commentType": "FUNCTION", "texts": [{"value": "Receptor for a ligand."}]},
            {"commentType": "PTM", "texts": [{"value": "Glycosylated."}]},
        ]
        p = normalize(raw)
        self.assertEqual(p["subunit_comments"][0]["text"], "Monomer. Forms a heterodimer with PARTNERX.")
        self.assertEqual(p["subunit_comments"][0]["eco"], ["ECO:0000269"])
        self.assertEqual(p["function_comments"][0]["text"], "Receptor for a ligand.")
        self.assertEqual(p["ptm_comments"][0]["text"], "Glycosylated.")
        unreviewed = normalize(dict(raw, entryType="UniProtKB unreviewed (TrEMBL)"))
        self.assertFalse(unreviewed["reviewed"])
        self.assertEqual(unreviewed["evidence"]["kind"], "prediction")


class MolecularProfileTests(unittest.TestCase):
    def test_profile_counts_and_positions(self):
        p = make_protein("P91001")
        p["features"] += [
            feature("disulfide", 22, 42, name="SS-internal"),
            feature("disulfide", 62, 90, name="SS-crossing"),
            feature("glycosylation_site", 32, 32, name="N-linked (GlcNAc...) asparagine"),
            feature("variant", 50, 50, name="R -> H", feature_id="VAR_9"),
            feature("mutagenesis", 55, 55, name="A->G: loss", feature_id="MUT_9"),
            feature("active_site", 88, 88, name="Proton acceptor"),
            feature("binding_site", 30, 30, name="Ligand"),
        ]
        candidates, _ = propose(p)
        primary = candidates[0]
        mp = primary["molecular_profile"]
        self.assertEqual(mp["length"], 66)
        self.assertEqual(mp["cysteine_count"], 3)
        self.assertTrue(mp["cysteine_odd"])
        self.assertEqual(mp["cysteine_positions"], [22, 42, 62])
        self.assertEqual([d["name"] for d in mp["disulfides_fully_contained"]], ["SS-internal"])
        self.assertEqual(mp["disulfides_partial"][0]["retained_position"], 62)
        self.assertEqual(mp["annotated_glycosylation_in_fragment"][0]["subtype"], "N-linked")
        self.assertEqual(len(mp["variants_overlapping"]), 1)
        self.assertEqual(mp["variants_overlapping"][0]["feature_id"], "VAR_9")
        self.assertEqual(len(mp["mutagenesis_records_overlapping"]), 1)
        self.assertEqual(mp["active_sites_in_fragment"], [])
        self.assertEqual(mp["active_sites_outside_fragment"][0]["start"], 88)
        self.assertEqual(mp["binding_sites_in_fragment"][0]["start"], 30)
        # Sequence-level sequon scan stays separate from annotated sites.
        self.assertEqual(mp["n_glyco_sequons"]["count"], 0)  # fixture N is always followed by P

    def test_sequon_scan_positions(self):
        self.assertEqual(n_glyco_sequons("ANATNPST", 10), [11])  # N-P does not count
        self.assertEqual(n_glyco_sequons("NQS", 1), [1])


class WarningLevelTests(unittest.TestCase):
    def test_free_thiol_odd_cysteine_is_warning_not_class_change(self):
        temp, result, records, candidates = screen_one(make_protein("P92001"))
        rec = records[0]
        self.assertEqual(rec["screening_recommendation"], "standard_candidate")
        primary = next(c for c in candidates if c.get("is_primary"))
        self.assertIn("free_thiol_odd_cysteine", codes(primary))
        self.assertEqual(primary["risks"] and "review",
                         next(r for r in primary["risks"] if r["code"] == "free_thiol_odd_cysteine")["severity"])
        self.assertNotIn("free_thiol_odd_cysteine", rec["reason_codes"])
        temp.cleanup()

    def test_even_cysteine_no_flag(self):
        p = make_protein("P92002", lambda p: p.update(sequence=p["sequence"][:21] + "S" + p["sequence"][22:]))
        candidates, _ = propose(p)
        self.assertNotIn("free_thiol_odd_cysteine", codes(candidates[0]))

    def test_dense_glycosylation_warning_not_class_change(self):
        def glyco_rich(p):
            # ECD (11-76) becomes sequon-dense: 20 sequons in 66 aa.
            p["sequence"] = p["sequence"][:10] + "NAT" * 20 + p["sequence"][70:]
        temp, result, records, candidates = screen_one(make_protein("P92003", glyco_rich))
        rec = records[0]
        primary = next(c for c in candidates if c.get("is_primary"))
        self.assertEqual(rec["screening_recommendation"], "standard_candidate")
        self.assertIn("dense_glycosylation", codes(primary))
        self.assertEqual(primary["molecular_profile"]["n_glyco_sequons"]["count"], 20)
        temp.cleanup()

    def test_dense_glycosylation_threshold_respected(self):
        # Few sequons in a long fragment: below the count heuristic -> no flag.
        p = make_protein("P92004", lambda p: p.update(sequence=p["sequence"][:10] + "NAT" + p["sequence"][13:]))
        candidates, _ = propose(p)
        self.assertNotIn("dense_glycosylation", codes(candidates[0]))
        self.assertIn("potential_n_glycosylation", codes(candidates[0]))


class DomainCutTests(unittest.TestCase):
    def test_domain_cut_by_boundary_is_conditional_with_identity(self):
        p = make_protein("P93001")
        p["features"].append(feature("domain", 70, 90, name="CROSS-ECD-TM-domain"))
        temp, result, records, candidates = screen_one(p)
        rec = records[0]
        self.assertEqual(rec["screening_recommendation"], "conditional_candidate")
        self.assertIn("domain_cut_by_boundary", rec["reason_codes"])
        self.assertIn("CROSS-ECD-TM-domain", rec["rationale"])
        primary = next(c for c in candidates if c.get("is_primary"))
        cut = primary["molecular_profile"]["domains_cut"]
        self.assertEqual(cut[0]["name"], "CROSS-ECD-TM-domain")
        self.assertEqual(cut[0]["cut_side"], "c_side")
        self.assertEqual(cut[0]["retained_interval"], [70, 76])
        temp.cleanup()

    def test_fully_contained_domains_not_cut(self):
        candidates, _ = propose(make_protein("P93002"))
        mp = candidates[0]["molecular_profile"]
        self.assertEqual({d["name"] for d in mp["domains_fully_contained"]},
                         {"synthetic-domain-A", "synthetic-domain-B"})
        self.assertEqual(mp["domains_cut"], [])


class BoundaryPrecisionTests(unittest.TestCase):
    def test_exact_adjacency_records_defining_feature(self):
        candidates, _ = propose(make_protein("P94001"))
        ba = candidates[0]["boundary_analysis"]
        self.assertEqual(ba["defining_feature"]["kind"], "extracellular")
        self.assertEqual(ba["n_side"]["status"], "exact")
        self.assertEqual(ba["n_side"]["adjacent_feature"]["kind"], "signal_peptide")
        self.assertEqual(ba["c_side"]["status"], "exact")
        self.assertEqual(ba["c_side"]["adjacent_feature"]["kind"], "transmembrane")
        self.assertEqual(ba["reason_codes"], [])

    def test_tm_adjacency_tolerance_used_is_informational(self):
        def shifted_tm(p):
            p["features"] = [f for f in p["features"] if f["kind"] != "transmembrane"]
            p["features"].append(feature("transmembrane", 79, 98))  # 2-residue unannotated gap
        temp, result, records, candidates = screen_one(make_protein("P94002", shifted_tm))
        rec = records[0]
        self.assertEqual(rec["screening_recommendation"], "standard_candidate")  # informational only
        self.assertIn("tm_adjacency_tolerance_used", rec["reason_codes"])
        primary = next(c for c in candidates if c.get("is_primary"))
        self.assertEqual(primary["boundary_analysis"]["c_side"]["status"], "gap_tolerated")
        temp.cleanup()

    def test_boundary_feature_conflict_is_conditional(self):
        def far_tm(p):
            p["features"] = [f for f in p["features"] if f["kind"] != "transmembrane"]
            p["features"].append(feature("transmembrane", 85, 100))  # gap of 8 > tolerance 5
        temp, result, records, candidates = screen_one(make_protein("P94003", far_tm))
        rec = records[0]
        self.assertEqual(rec["screening_recommendation"], "conditional_candidate")
        self.assertIn("boundary_feature_conflict", rec["reason_codes"])
        self.assertIn("c_side", rec["rationale"])
        primary = next(c for c in candidates if c.get("is_primary"))
        self.assertEqual(primary["boundary_analysis"]["c_side"]["status"], "gap_beyond_tolerance")
        temp.cleanup()

    def test_missing_adjacent_feature_recorded_not_flagged(self):
        def no_sp(p):
            p["features"] = [f for f in p["features"] if f["kind"] != "signal_peptide"]
        candidates, _ = propose(make_protein("P94004", no_sp))
        ba = candidates[0]["boundary_analysis"]
        self.assertEqual(ba["n_side"]["status"], "missing")
        self.assertNotIn("boundary_feature_conflict", codes(candidates[0]))
        # The existing SP-missing review flag still fires.
        self.assertIn("signal_peptide_annotation_missing", codes(candidates[0]))

    def test_domain_alternative_adjacency_not_forced(self):
        p = make_protein("P94005")
        p["candidate_policy"] = "full_ecd_and_domains"
        candidates, _ = propose(p)
        domain_cand = next(c for c in candidates if c["antigen_form_type"] == "domain_fragment")
        self.assertEqual(domain_cand["boundary_analysis"]["n_side"]["status"], "not_applicable")


class MultichainTests(unittest.TestCase):
    def test_heterooligomer_partner_is_conditional_with_named_evidence(self):
        def subunit(p):
            p["subunit_comments"] = [{"text": "Monomer. Forms a heterodimer with PARTNERX; stable at the surface.", "eco": []}]
        temp, result, records, candidates = screen_one(make_protein("P95001", subunit))
        rec = records[0]
        self.assertEqual(rec["screening_recommendation"], "conditional_candidate")
        self.assertIn("native_heteromer_context_requires_review", rec["reason_codes"])
        self.assertIn("heterodimer with PARTNERX", rec["rationale"])
        self.assertEqual(rec["evidence"]["multichain_partners"],
                         ["Forms a heterodimer with PARTNERX; stable at the surface."])
        temp.cleanup()

    def test_homodimer_and_plain_interaction_not_flagged(self):
        self.assertEqual(detect_multichain_partners(
            [{"text": "Homodimer; disulfide-linked. Interacts with SH3BP3."}]), [])
        self.assertEqual(detect_multichain_partners(
            [{"text": "Part of a complex composed of A, B and C."}]), [])  # documented conservative limit

    def test_no_subunit_comment_no_flag(self):
        candidates, _ = propose(make_protein("P95003"))
        self.assertNotIn("native_heteromer_context_requires_review", codes(candidates[0]))
        self.assertEqual(candidates[0]["multichain_partners"], [])


class MultiPassLoopTests(unittest.TestCase):
    def _multi(self, pid, with_domain=True):
        def build(p):
            p.update(topology="multi_pass")
            feats = [feature("extracellular", 10, 20), feature("transmembrane", 21, 39),
                     feature("extracellular", 40, 50), feature("transmembrane", 51, 60)]
            if with_domain:
                feats.append(feature("domain", 10, 20, name="ext-loop-1-domain"))
            p["features"] = feats
        return make_protein(pid, build)

    def test_loop_enumeration_explicit_with_alternate_link(self):
        temp, result, records, candidates = screen_one(self._multi("P96001"))
        rec = records[0]
        loops = rec["multipass_loops"]
        self.assertEqual(len(loops), 2)
        self.assertEqual([(l["start"], l["end"], l["length"]) for l in loops],
                         [(10, 20, 11), (40, 50, 11)])
        self.assertTrue(loops[0]["eligible_domain_alternate"])
        self.assertEqual(loops[0]["domains_fully_contained"][0]["name"], "ext-loop-1-domain")
        self.assertTrue(loops[0]["alternate_candidate_ids"])  # linked to the proposed alternate
        self.assertFalse(loops[1]["eligible_domain_alternate"])
        self.assertEqual(rec["screening_recommendation"], "conditional_candidate")
        temp.cleanup()

    def test_loop_enumeration_recorded_even_without_alternate(self):
        temp, result, records, candidates = screen_one(self._multi("P96002", with_domain=False))
        rec = records[0]
        self.assertEqual(rec["screening_recommendation"], "no_standard_route")
        loops = rec["multipass_loops"]
        self.assertEqual(len(loops), 2)
        self.assertFalse(any(l["eligible_domain_alternate"] for l in loops))
        self.assertTrue(all(l["alternate_candidate_ids"] == [] for l in loops))
        self.assertTrue(any("逐胞外环枚举" in m for m in rec["missing_info"]))
        temp.cleanup()

    def test_enumerate_external_loops_skips_fuzzy_loop_honestly(self):
        p = self._multi("P96003")
        fuzzy = feature("extracellular", 70, None, name="fuzzy-loop")
        fuzzy["boundary_status"] = "fuzzy"
        loops = enumerate_external_loops(p, p["features"] + [fuzzy], len(p["sequence"]))
        self.assertEqual(loops[-1]["length"], None)
        self.assertIn("deferred", loops[-1]["note"])

    def test_multipass_nonstandard_residue_blocks_gracefully_not_crash(self):
        """v0.4.1 regression (found by the 7,783-entry full run, Q9C0D9/SELENOI):
        a multi-pass protein carrying selenocysteine (U) must be classified
        with the reference_sequence_invalid block from propose(); the loop
        enumeration must not re-raise on the invalid sequence."""
        p = self._multi("P96004")
        p["sequence"] = p["sequence"][:30] + "U" + p["sequence"][31:]
        temp, result, records, candidates = screen_one(p)
        rec = records[0]
        self.assertEqual(rec["processing_status"], "ok")
        self.assertIn("reference_sequence_invalid", rec["reason_codes"])
        self.assertNotIn("multipass_loops", rec)  # skipped honestly, not faked
        # Invalid reference sequence means the entry cannot be evaluated, which
        # the decision table maps to insufficient_evidence (not a route verdict).
        self.assertEqual(rec["screening_recommendation"], "insufficient_evidence")
        temp.cleanup()


class ScopedPredictionTests(unittest.TestCase):
    def test_predicted_informational_feature_no_longer_raises_boundary_flag(self):
        def predicted_glyco(p):
            g = feature("glycosylation_site", 32, 32, name="N-linked (GlcNAc...) asparagine")
            g["evidence"] = dict(g["evidence"], kind="prediction")
            p["features"].append(g)
        candidates, _ = propose(make_protein("P97001", predicted_glyco))
        self.assertNotIn("prediction_requires_annotation_review", codes(candidates[0]))

    def test_predicted_boundary_feature_still_raises(self):
        def predicted_ecd(p):
            for f in p["features"]:
                if f["kind"] == "extracellular":
                    f["evidence"] = dict(f["evidence"], kind="prediction")
        candidates, _ = propose(make_protein("P97002", predicted_ecd))
        self.assertIn("prediction_requires_annotation_review", codes(candidates[0]))

    def test_invalid_informational_feature_deferred_not_blocking(self):
        def fuzzy_variant(p):
            v = feature("variant", 50, None, name="fuzzy variant")
            v["boundary_status"] = "fuzzy"
            p["features"].append(v)
        temp, result, records, candidates = screen_one(make_protein("P97003", fuzzy_variant))
        rec = records[0]
        self.assertEqual(rec["screening_recommendation"], "standard_candidate")
        self.assertEqual(rec["deferred_annotations"][0]["kind"], "variant")
        temp.cleanup()


class OutputContractTests(unittest.TestCase):
    def test_screening_outputs_carry_profile_and_boundary_fields(self):
        temp, result, records, candidates = screen_one(make_protein("P98001"), pilot=True)
        rec = records[0]
        primary = next(c for c in candidates if c.get("is_primary"))
        self.assertIn("molecular_profile", primary)
        self.assertIn("boundary_analysis", primary)
        self.assertIn("molecular_profile", rec["evidence"])
        self.assertIn("boundary_analysis", rec["evidence"])
        self.assertIn("reviewed", rec["evidence"]["reference"])
        self.assertIn("annotation_score", rec["evidence"]["reference"])
        summary = result["summary"]
        self.assertIn("reason_code_tallies", summary)
        tallies = summary["reason_code_tallies"]["primary_candidate_review_warnings"]
        self.assertEqual(tallies["denominator"], 1)
        self.assertEqual(tallies["counts"].get("free_thiol_odd_cysteine"), 1)
        pi = (Path(result["run_dir"]) / "PI_SUMMARY.md").read_text()
        self.assertIn("checks applied", pi)
        self.assertIn("已评估", pi)
        self.assertIn("未评估", pi)
        tsv_text = (Path(result["run_dir"]) / "protein_screening.tsv").read_text()
        self.assertIn("primary_molecular_profile", tsv_text.splitlines()[0])
        self.assertIn("multipass_loops", tsv_text.splitlines()[0])
        temp.cleanup()


if __name__ == "__main__":
    unittest.main()
