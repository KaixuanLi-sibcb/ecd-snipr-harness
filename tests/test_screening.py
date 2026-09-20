import copy
import json
import tempfile
import unittest
from pathlib import Path
from test_design import ROOT, approve, feature, project, protein, scaffold
from ecd_snipr.common import file_hash, read_json, write_json
from ecd_snipr.acquisition import (build_set, disposition, map_gene_name, read_target_list,
                                   search_accessions, select_analysis_reference)
from ecd_snipr.design import propose
from ecd_snipr.harness import run as harness_run, verify_bundle
from ecd_snipr.screening import run_screening


def make_protein(pid, mutate=None):
    p = copy.deepcopy(project()["proteins"][0])
    p["protein_id"] = pid
    p["accession"] = pid
    p["candidate_policy"] = "full_ecd_only"
    p["candidate_regions"] = []
    if mutate:
        mutate(p)
    select_analysis_reference(p)
    return p


def make_set(setdir, proteins, extra_entries=(), failed_rows=()):
    """Write a minimal set_definition.json plus per-protein files (offline)."""
    setdir = Path(setdir)
    (setdir / "proteins").mkdir(parents=True, exist_ok=True)
    entries = []
    for i, p in enumerate(proteins):
        write_json(setdir / "proteins" / (p["accession"] + ".json"), p)
        try:
            disp = disposition(p)
        except Exception:
            disp = {"scope_status": "unevaluated", "note": "disposition failed on malformed record"}
        entries.append({"entry_id": f"entry-{i + 1:05d}", "input": {"row": i + 1, "value": p["accession"], "input_type": "accession"},
                        "identity_status": "resolved", "accession": p["accession"], "gene": p.get("gene", ""),
                        "mapping": None, "fetch": None,
                        "protein_file": "proteins/" + p["accession"] + ".json",
                        "protein_sha256": file_hash(setdir / "proteins" / (p["accession"] + ".json")),
                        "reference_selection": p["reference_selection"],
                        "processing_status": "ok", "processing_error": "", "disposition": disp})
    entries.extend(extra_entries)
    n = len(entries)
    definition = {"set_id": "set-synthetic", "created_at": "2026-01-01T00:00:00Z", "generator_version": "test",
                  "database": {"name": "UniProtKB", "access": "test", "release": "synthetic", "query_record": None},
                  "query": None, "inclusion_rules": ["synthetic test set"],
                  "counts": {"input_rows": n, "resolved": len(proteins),
                             "unresolved": n - len(proteins), "technical_failures": len(failed_rows),
                             "duplicate_input_rows": sum(bool(e.get("duplicate_of")) for e in entries)},
                  "completeness": "partial" if failed_rows else "complete", "entries": entries}
    write_json(setdir / "set_definition.json", definition)
    return definition


def to_secreted(p):
    """A realistic secreted record: no TM/cytoplasmic/topological-domain annotation."""
    p.update(topology="secreted", location="secreted")
    p["features"] = [f for f in p["features"] if f["kind"] in {"signal_peptide", "domain", "epitope", "disulfide"}]
    p["features"].append(feature("chain", 11, 95))


class Raw:
    """Minimal UniProt-shaped payload for offline acquisition tests."""
    @staticmethod
    def entry(accession="P00001", gene="SYNTH", description="Extracellular", tm=(71, 90), ext=(11, 70), version=3):
        return {"primaryAccession": accession, "entryType": "UniProtKB reviewed (Swiss-Prot)",
                "sequence": {"value": "ACDEFGHIKLMNPQRSTVWY" * 5}, "organism": {"taxonId": 9606},
                "entryAudit": {"entryVersion": version},
                "genes": [{"geneName": {"value": gene}}],
                "comments": [{"commentType": "SUBCELLULAR LOCATION",
                              "subcellularLocations": [{"location": {"value": "Cell membrane"}}]}],
                "features": [
                    {"type": "Topological domain", "description": description,
                     "location": {"start": {"value": ext[0]}, "end": {"value": ext[1]}}},
                    {"type": "Transmembrane",
                     "location": {"start": {"value": tm[0]}, "end": {"value": tm[1]}}}]}


class Response:
    def __init__(self, payload, headers=None):
        self._payload = json.dumps(payload).encode()
        self.headers = headers or {}

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def read(self):
        return self._payload


class AcquisitionTests(unittest.TestCase):
    def test_target_list_preserves_every_row(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "list.tsv"
            path.write_text("accession\tgene\nP15391\tCD19\n\tLONESOME\nQ15116\tPDCD1\nBAD-ACC-XX\tFOO\n\nP18627\t\n", encoding="utf-8")
            rows = read_target_list(path)
            self.assertEqual(len(rows), 6)  # blank row preserved
            self.assertEqual(rows[0]["input_type"], "accession")
            self.assertEqual(rows[1]["input_type"], "gene_name")
            self.assertEqual(rows[3]["input_type"], "invalid_accession")
            self.assertEqual(rows[4]["input_type"], "blank")
            self.assertEqual(rows[5]["input_type"], "accession")
            self.assertEqual(rows[1]["row"], 3)

    def test_gene_mapping_records_evidence_and_ambiguity(self):
        two = {"results": [Raw.entry("P00001", "DUP"), Raw.entry("P00002", "DUP")]}
        one = {"results": [Raw.entry("P00003", "SOLO")]}
        none = {"results": []}
        with tempfile.TemporaryDirectory() as temp:
            opener = lambda *a, **k: Response(two, {"X-UniProt-Release": "2026_01"})
            mapping = map_gene_name("dup", temp, opener=opener, sleeper=lambda _: None)
            self.assertEqual(mapping["status"], "ambiguous")
            self.assertEqual({c["accession"] for c in mapping["candidates"]}, {"P00001", "P00002"})
            self.assertTrue(mapping["queries"])
            cached = map_gene_name("dup", temp, offline=True)
            self.assertEqual(cached["status"], "ambiguous")  # cache reuse offline
            mapping = map_gene_name("solo", temp, opener=lambda *a, **k: Response(one), sleeper=lambda _: None)
            self.assertEqual((mapping["status"], mapping["accession"]), ("resolved", "P00003"))
            mapping = map_gene_name("ghost", temp, opener=lambda *a, **k: Response(none), sleeper=lambda _: None)
            self.assertEqual(mapping["status"], "unresolved")
            self.assertEqual(len(mapping["queries"]), 2)  # reviewed then broader query recorded
            mapping = map_gene_name("neverfetched", temp, offline=True)
            self.assertEqual(mapping["status"], "missing")

    def test_search_accessions_records_release_and_limit(self):
        page = {"results": [Raw.entry(f"P{i:05d}", f"G{i}") for i in range(3)]}
        headers = {"X-UniProt-Release": "2026_02", "X-Total-Results": "10"}
        result = search_accessions("(organism_id:9606)", limit=3, opener=lambda *a, **k: Response(page, headers), sleeper=lambda _: None)
        self.assertEqual(result["status"], "fetched")
        self.assertEqual(result["uniprot_release"], "2026_02")
        self.assertEqual(result["completeness"], "truncated_by_limit")
        self.assertEqual(len(result["entries"]), 3)

    def test_build_set_single_failure_isolated_and_partial(self):
        raw_ok = Raw.entry("P00001", "OKG")
        def opener(request, **kwargs):
            url = request.full_url if hasattr(request, "full_url") else str(request)
            if "P00002" in url:
                raise OSError("simulated fetch failure")
            return Response(raw_ok if "P00001" in url else Raw.entry("P00003", "G3"))
        with tempfile.TemporaryDirectory() as temp:
            rows = [{"row": i + 1, "value": a, "input_type": "accession", "accession": a, "gene_name": ""}
                    for i, a in enumerate(("P00001", "P00002", "P00003"))]
            definition = build_set(rows=rows, cache=Path(temp) / "cache", outdir=Path(temp) / "set",
                                   opener=opener, sleeper=lambda _: None)
            self.assertEqual(definition["counts"]["resolved"], 2)
            self.assertEqual(definition["counts"]["technical_failures"], 1)
            self.assertEqual(definition["completeness"], "partial")
            failed = next(e for e in definition["entries"] if e["accession"] == "P00002")
            self.assertEqual(failed["processing_status"], "failed")
            self.assertIn("simulated fetch failure", failed["processing_error"])

    def test_build_set_duplicate_rows_linked_not_dropped(self):
        raw = Raw.entry("P00001", "DUPG")
        with tempfile.TemporaryDirectory() as temp:
            rows = [{"row": i + 1, "value": "P00001", "input_type": "accession", "accession": "P00001", "gene_name": ""}
                    for i in range(2)]
            definition = build_set(rows=rows, cache=Path(temp) / "c", outdir=Path(temp) / "s",
                                   opener=lambda *a, **k: Response(raw), sleeper=lambda _: None)
            self.assertEqual(len(definition["entries"]), 2)
            self.assertEqual(definition["entries"][1].get("duplicate_of"), "entry-00001")
            self.assertEqual(definition["counts"]["duplicate_input_rows"], 1)

    def test_canonical_selection_is_not_experimental_confirmation(self):
        from ecd_snipr.uniprot import normalize
        fresh = normalize(Raw.entry("P00001"))  # plain canonical fetch: isoform unset, ambiguity flagged
        self.assertTrue(fresh["isoform_ambiguous"])
        p = select_analysis_reference(fresh)
        self.assertEqual(p["reference_selection"]["basis"], "database_canonical")
        self.assertEqual(p["reference_selection"]["experimental_isoform_confirmation"], "not_performed")
        self.assertFalse(p["isoform_ambiguous"])
        self.assertEqual(p["isoform"], "P00001:canonical")
        # Screening proceeds, but the missing confirmation is explicit.
        candidates, flags = propose(p)
        self.assertNotIn("isoform_unresolved", {f["code"] for f in flags})
        # Genuine ambiguity survives: caller keeps isoform_ambiguous and screening blocks.
        q = dict(p)
        q["isoform"], q["isoform_ambiguous"] = "", True
        del q["reference_selection"]
        _, flags = propose(q)
        self.assertIn("isoform_unresolved", {f["code"] for f in flags})

    def test_disposition_distinguishes_membership_surface_scope(self):
        secreted = make_protein("P00010", to_secreted)
        d = disposition(secreted)
        self.assertEqual((d["extension_set"], d["scope_status"]), ("secreted_extension", "in_scope"))
        self.assertEqual(d["natural_cell_surface_target"], "no")
        organelle = make_protein("P00011", lambda p: (p.update(location="other_membrane"),
                                                      p["features"].__setitem__(1, feature("lumenal", 11, 76))))
        d = disposition(organelle)
        self.assertEqual(d["scope_status"], "out_of_scope")
        self.assertIn("organelle_membrane_not_default_scope", d["scope_reason_codes"])
        multi = make_protein("P00012", lambda p: (p.update(topology="multi_pass"),
                                                  p.update(features=[feature("extracellular", 10, 20), feature("transmembrane", 21, 39),
                                                                     feature("extracellular", 40, 50), feature("transmembrane", 51, 60)])))
        d = disposition(multi)
        self.assertEqual(d["scope_status"], "in_scope")
        self.assertTrue(any("multi_pass" in n for n in d["disposition_notes"]))

    def test_secreted_location_precedence_over_generic_membrane(self):
        # Regression (v0.3.0 lab cross-check): UniProt entries annotated
        # "Secreted" plus a generic organelle membrane (REN: Secreted + Membrane;
        # TFPI: Secreted + Microsome membrane) were misread as other_membrane
        # and ruled out_of_scope. Secreted must win over generic membrane, while
        # plasma membrane still wins over both.
        from ecd_snipr.uniprot import normalize
        raw = Raw.entry()
        del raw["features"]  # no TM: secreted path
        raw["features"] = [{"type": "Signal", "location": {"start": {"value": 1}, "end": {"value": 10}}},
                           {"type": "Chain", "description": "Mature", "location": {"start": {"value": 11}, "end": {"value": 100}}}]
        raw["comments"] = [{"commentType": "SUBCELLULAR LOCATION",
                            "subcellularLocations": [{"location": {"value": "Secreted"}},
                                                     {"location": {"value": "Membrane"}}]}]
        p = normalize(raw)
        self.assertEqual(p["location"], "secreted")
        self.assertEqual(p["topology"], "secreted")
        d = disposition(p)
        self.assertEqual((d["scope_status"], d["extension_set"]), ("in_scope", "secreted_extension"))
        self.assertEqual(p["raw_locations"], ["Secreted", "Membrane"])  # raw values preserved

        raw["comments"] = [{"commentType": "SUBCELLULAR LOCATION",
                            "subcellularLocations": [{"location": {"value": "Secreted"}},
                                                     {"location": {"value": "Cell membrane"}}]}]
        self.assertEqual(normalize(raw)["location"], "plasma_membrane")

        raw["comments"] = [{"commentType": "SUBCELLULAR LOCATION",
                            "subcellularLocations": [{"location": {"value": "Microsome membrane"}}]}]
        self.assertEqual(normalize(raw)["location"], "other_membrane")


class ScreeningDecisionTests(unittest.TestCase):
    def _screen(self, proteins, extra_entries=(), failed_rows=(), **kwargs):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        setdir = Path(temp.name) / "set"
        make_set(setdir, proteins, extra_entries, failed_rows)
        result = run_screening(setdir, Path(temp.name) / "screen", **kwargs)
        records = read_json(Path(result["run_dir"]) / "screening.json")
        return result, records

    def test_screening_completes_without_scaffold_review_experiments(self):
        result, records = self._screen([make_protein("P00001")], pilot=True)
        rec = records[0]
        self.assertEqual(rec["screening_recommendation"], "standard_candidate")
        candidates = read_json(Path(result["run_dir"]) / "candidates.json")
        self.assertTrue(candidates)
        self.assertTrue(all(not c["fusion_sequence"] for c in candidates))
        self.assertTrue(all(c["functional_status"] == "not_experimentally_validated" for c in candidates))
        fasta = (Path(result["run_dir"]) / "candidate_fragments.fasta").read_text()
        self.assertIn("UNVALIDATED", fasta)
        self.assertTrue(candidates[0]["junction_notes"])

    def test_missing_optional_risk_info_does_not_downgrade(self):
        _, records = self._screen([make_protein("P00001")])
        rec = records[0]
        self.assertEqual(rec["screening_recommendation"], "standard_candidate")
        self.assertTrue(any("上下文风险记录缺失" in m for m in rec["missing_info"]))
        self.assertTrue(any("isoform" in m for m in rec["missing_info"]))
        self.assertIn("positive_sequence_topology_boundary_evidence", rec["reason_codes"])

    def test_core_sequence_errors_still_blocked(self):
        bad = make_protein("P00020")
        bad["sequence"] = bad["sequence"] + "X"
        _, records = self._screen([bad])
        self.assertEqual(records[0]["screening_recommendation"], "insufficient_evidence")
        self.assertIn("reference_sequence_invalid", records[0]["reason_codes"])
        # A proposed region retaining native TM is blocked; the valid baseline still wins.
        p = make_protein("P00021")
        p["candidate_regions"] = [dict(feature("region", 11, 82), rationale="Deliberately invalid test region")]
        result, records = self._screen([p])
        candidates = read_json(Path(result["run_dir"]) / "candidates.json")
        blocked = [c for c in candidates if c["design_status"] == "blocked"]
        self.assertTrue(any("retained_transmembrane" in {r["code"] for r in c["risks"]} for c in blocked))
        self.assertEqual(records[0]["screening_recommendation"], "standard_candidate")

    def test_type_ii_records_orientation_issue(self):
        p = make_protein("P00030", lambda p: (p.update(topology="type_ii"),
                                              p.update(features=[feature("cytoplasmic", 1, 10), feature("transmembrane", 11, 20),
                                                                 feature("extracellular", 21, 95), feature("domain", 30, 80)])))
        _, records = self._screen([p])
        rec = records[0]
        self.assertEqual(rec["screening_recommendation"], "conditional_candidate")
        self.assertIn("type_ii_attachment_orientation_change", rec["reason_codes"])

    def test_gpi_mature_boundary_checked(self):
        def gpi(p):
            p.update(topology="gpi")
            p["features"] = [feature("signal_peptide", 1, 10), feature("extracellular", 11, 80),
                             feature("gpi_attachment_site", 80, 80), feature("gpi_signal", 81, 95),
                             feature("chain", 11, 80)]
        _, records = self._screen([make_protein("P00031", gpi)])
        rec = records[0]
        self.assertEqual(rec["screening_recommendation"], "conditional_candidate")
        self.assertIn("gpi_anchor_replaced", rec["reason_codes"])

    def test_multipass_domain_alternates_vs_no_route(self):
        def multi(p):
            p.update(topology="multi_pass")
            p["features"] = [feature("extracellular", 10, 20), feature("transmembrane", 21, 39),
                             feature("extracellular", 40, 50), feature("transmembrane", 51, 60),
                             feature("domain", 10, 20, name="ext-loop-1")]
        _, records = self._screen([make_protein("P00032", multi)])
        rec = records[0]
        self.assertEqual(rec["screening_recommendation"], "conditional_candidate")
        self.assertIn("discontinuous_or_multipass", rec["reason_codes"])
        self.assertTrue(rec["alternates"] or rec["primary_candidate_id"])
        def multi_bare(p):
            p.update(topology="multi_pass")
            p["features"] = [feature("extracellular", 10, 20), feature("transmembrane", 21, 39),
                             feature("extracellular", 40, 50), feature("transmembrane", 51, 60)]
        _, records = self._screen([make_protein("P00033", multi_bare)])
        rec = records[0]
        self.assertEqual(rec["screening_recommendation"], "no_standard_route")
        self.assertIn("路线限制", rec["rationale"])
        self.assertNotIn("永远", rec["rationale"].replace("不表示该蛋白永远不可用", ""))

    def test_secreted_extension_counted_separately(self):
        p = make_protein("P00034", to_secreted)
        result, records = self._screen([p])
        self.assertEqual(records[0]["screening_recommendation"], "conditional_candidate")
        counts = result["summary"]["counts"]
        self.assertEqual(counts["secreted_extension_references"], 1)
        self.assertEqual(counts["core_references"], 0)
        self.assertEqual(counts["screening_classes_core_only"]["conditional_candidate"], 0)

    def test_intracellular_out_of_scope_not_a_class(self):
        p = make_protein("P00035", lambda p: p.update(topology="intracellular", location="unknown", features=[]))
        result, records = self._screen([p])
        rec = records[0]
        self.assertIsNone(rec["screening_recommendation"])
        self.assertEqual(rec["scope_status"], "out_of_scope")
        self.assertEqual(result["summary"]["counts"]["out_of_scope"], 1)
        self.assertEqual(result["summary"]["counts"]["evaluated_references"], 0)

    def test_identity_ambiguity_preserved(self):
        ambiguous = {"entry_id": "entry-00002", "input": {"row": 2, "value": "DUPGENE", "input_type": "gene_name"},
                     "identity_status": "ambiguous", "accession": "", "gene": "DUPGENE",
                     "mapping": {"mapping_status": "ambiguous", "candidates": [{"accession": "P00001"}, {"accession": "P00002"}]},
                     "fetch": None, "protein_file": "", "protein_sha256": "", "reference_selection": None,
                     "processing_status": "ok", "processing_error": "",
                     "disposition": {"scope_status": "unevaluated", "note": "identity not resolved; ambiguity preserved"}}
        _, records = self._screen([make_protein("P00001")], extra_entries=[ambiguous])
        rec = next(r for r in records if r["entry_id"] == "entry-00002")
        self.assertEqual(rec["screening_recommendation"], "insufficient_evidence")
        self.assertIn("identity_ambiguous", rec["reason_codes"])
        self.assertEqual(len(rec["evidence"]["mapping"]["candidates"]), 2)

    def test_lab_rules_absent_marker_and_provenance(self):
        _, records = self._screen([make_protein("P00001")])
        self.assertEqual(records[0]["lab_rules"]["status"], "尚未纳入")
        with tempfile.TemporaryDirectory() as temp:
            rules = {"rules": [{"rule_id": "LAB-1", "match": {"accession": "P00001"}, "effect": "downgrade_to_conditional",
                                "reason_code": "lab_expression_difficult", "rationale": "实验室记录该靶点表达困难",
                                "evidence": {"kind": "local_experiment", "source": "lab-notebook-7", "version": "2026-01"}}]}
            path = Path(temp) / "rules.json"
            write_json(path, rules)
            _, records = self._screen([make_protein("P00001")], lab_rules_path=path)
            rec = records[0]
            self.assertEqual(rec["screening_recommendation"], "conditional_candidate")
            self.assertEqual(rec["lab_rules"]["status"], "applied")
            self.assertEqual(rec["lab_rules"]["applied"][0]["evidence"]["source"], "lab-notebook-7")
            # Rules never upgrade: a conditional target stays conditional without rules.
            def multi(p):
                p.update(topology="multi_pass")
                p["features"] = [feature("extracellular", 10, 20), feature("transmembrane", 21, 39),
                                 feature("extracellular", 40, 50), feature("transmembrane", 51, 60),
                                 feature("domain", 10, 20, name="loop")]
            upgrade = {"rules": [{"rule_id": "LAB-2", "match": {"topology": "multi_pass"}, "effect": "annotate",
                                  "reason_code": "note_only", "rationale": "注释，不升级",
                                  "evidence": {"kind": "local_experiment", "source": "nb", "version": "1"}}]}
            write_json(path, upgrade)
            _, records = self._screen([make_protein("P00036", multi)], lab_rules_path=path)
            self.assertEqual(records[0]["screening_recommendation"], "conditional_candidate")
            self.assertIn("note_only", records[0]["reason_codes"])
            bad = {"rules": [{"rule_id": "X", "match": {"gene": "G"}, "effect": "downgrade_to_conditional",
                              "rationale": "no provenance", "evidence": {"kind": "prediction"}}]}
            write_json(path, bad)
            with self.assertRaises(ValueError):
                self._screen([make_protein("P00037")], lab_rules_path=path)

    def test_primary_selection_and_alternate_reasons(self):
        p = make_protein("P00040")
        p["candidate_policy"] = "full_ecd_and_domains"
        result, records = self._screen([p])
        rec = records[0]
        candidates = read_json(Path(result["run_dir"]) / "candidates.json")
        primary = next(c for c in candidates if c["candidate_id"] == rec["primary_candidate_id"])
        self.assertEqual(primary["antigen_form_type"], "full_ecd")
        self.assertTrue(rec["alternates"])
        self.assertTrue(all("reason_not_primary" in a for a in rec["alternates"]))
        self.assertTrue(any(not c["is_primary"] for c in candidates))

    def test_multiple_candidates_each_evaluated(self):
        def multichain(p):
            to_secreted(p)
            p["features"] = [f for f in p["features"] if f["kind"] != "chain"] + [feature("chain", 11, 50), feature("chain", 51, 95)]
        result, records = self._screen([make_protein("P00041", multichain)])
        rec = records[0]
        self.assertGreaterEqual(rec["candidates_evaluated"], 2)
        self.assertEqual(rec["screening_recommendation"], "conditional_candidate")
        self.assertIn("multiple_processed_chains_dependency_unknown", rec["reason_codes"])


    def test_fuzzy_secondary_annotation_deferred_not_blocking(self):
        # LAG3-like: exact ECD/chain plus a fuzzy shed-form chain annotation.
        p = make_protein("P00042")
        fuzzy_chain = feature("chain", 23, None, name="shed-form")
        fuzzy_chain["boundary_status"] = "fuzzy"
        p["features"].append(fuzzy_chain)
        result, records = self._screen([p])
        rec = records[0]
        self.assertEqual(rec["screening_recommendation"], "standard_candidate")
        self.assertTrue(rec.get("deferred_annotations"))
        self.assertEqual(rec["deferred_annotations"][0]["kind"], "chain")
        self.assertIn("secondary_feature_annotation_deferred", rec["reason_codes"])
        # Genuine essential coordinate problems still block: sole fuzzy SP.
        q = make_protein("P00043")
        q["features"] = [f for f in q["features"] if f["kind"] != "signal_peptide"]
        fuzzy_sp = feature("signal_peptide", 1, None)
        fuzzy_sp["boundary_status"] = "fuzzy"
        q["features"].insert(0, fuzzy_sp)
        _, records = self._screen([q])
        self.assertEqual(records[0]["screening_recommendation"], "insufficient_evidence")
        self.assertIn("feature_annotation_invalid", records[0]["reason_codes"])
        # Secreted route: a fuzzy sole chain is essential and still blocks.
        s = make_protein("P00044", to_secreted)
        s["features"] = [f for f in s["features"] if f["kind"] != "chain"]
        fuzzy_only_chain = feature("chain", 11, None)
        fuzzy_only_chain["boundary_status"] = "fuzzy"
        s["features"].append(fuzzy_only_chain)
        _, records = self._screen([s])
        self.assertEqual(records[0]["screening_recommendation"], "insufficient_evidence")

    def test_prediction_boundaries_reviewed_vs_unreviewed(self):
        reviewed = make_protein("P00045")
        reviewed["evidence"] = {"kind": "curated_annotation", "source": "synthetic://software-tests", "version": "1"}
        reviewed["features"][1]["evidence"]["kind"] = "prediction"  # ECD boundary predicted
        _, records = self._screen([reviewed])
        self.assertEqual(records[0]["screening_recommendation"], "standard_candidate")  # reviewed entry: risk stays review-level
        unreviewed = make_protein("P00046")
        unreviewed["features"][1]["evidence"]["kind"] = "prediction"
        unreviewed["evidence"] = {"kind": "prediction", "source": "synthetic://software-tests", "version": "1"}
        _, records = self._screen([unreviewed])
        self.assertEqual(records[0]["screening_recommendation"], "standard_candidate")
        self.assertNotIn("prediction_requires_annotation_review", records[0]["reason_codes"])


class RunnerTests(unittest.TestCase):
    def test_bundle_outputs_resume_and_verify(self):
        with tempfile.TemporaryDirectory() as temp:
            setdir = Path(temp) / "set"
            make_set(setdir, [make_protein("P00001")])
            result = run_screening(setdir, Path(temp) / "screen", pilot=True)
            run_dir = Path(result["run_dir"])
            for name in ("protein_screening.tsv", "candidate_plan.tsv", "candidates.json",
                         "candidate_fragments.fasta", "summary.json", "PI_SUMMARY.md",
                         "screening_overview.svg", "screening_overview_data.json",
                         "screening.json", "set_definition.json", "manifest.json"):
                self.assertTrue((run_dir / name).is_file(), name)
            self.assertTrue(verify_bundle(run_dir))
            self.assertEqual(run_screening(setdir, Path(temp) / "screen", pilot=True, resume=True)["execution"],
                             "verified_cache_hit")

    def test_partial_run_not_masquerading(self):
        failed_entry = {"entry_id": "entry-00002", "input": {"row": 2, "value": "P99999", "input_type": "accession"},
                        "identity_status": "unresolved", "accession": "P99999", "gene": "",
                        "mapping": None, "fetch": {"status": "error", "reason": "simulated"},
                        "protein_file": "", "protein_sha256": "", "reference_selection": None,
                        "processing_status": "failed", "processing_error": "OSError: simulated",
                        "disposition": {"scope_status": "unevaluated", "note": "technical failure during acquisition"}}
        with tempfile.TemporaryDirectory() as temp:
            setdir = Path(temp) / "set"
            make_set(setdir, [make_protein("P00001")], extra_entries=[failed_entry], failed_rows=[failed_entry])
            result = run_screening(setdir, Path(temp) / "screen")
            summary = result["summary"]
            self.assertEqual(summary["completeness"]["state"], "partial")
            self.assertEqual(summary["counts"]["technical_failures"], 1)
            self.assertTrue(verify_bundle(result["run_dir"]))  # integrity ok, completeness partial
            records = read_json(Path(result["run_dir"]) / "screening.json")
            failed = next(r for r in records if r["entry_id"] == "entry-00002")
            self.assertEqual(failed["processing_status"], "technical_failure")
            self.assertIsNone(failed["screening_recommendation"])
            self.assertIn("partial", (Path(result["run_dir"]) / "PI_SUMMARY.md").read_text().split("集合完整性：")[1])

    def test_screening_crash_isolated_per_item(self):
        broken = make_protein("P00050")
        broken["features"] = "not-a-list"  # propose() will raise; batch must continue
        with tempfile.TemporaryDirectory() as temp:
            setdir = Path(temp) / "set"
            make_set(setdir, [broken, make_protein("P00051")])
            result = run_screening(setdir, Path(temp) / "screen")
            records = read_json(Path(result["run_dir"]) / "screening.json")
            crashed = next(r for r in records if r["accession"] == "P00050")
            self.assertEqual(crashed["processing_status"], "technical_failure")
            ok = next(r for r in records if r["accession"] == "P00051")
            self.assertEqual(ok["screening_recommendation"], "standard_candidate")

    def test_summary_ratios_have_denominators(self):
        with tempfile.TemporaryDirectory() as temp:
            setdir = Path(temp) / "set"
            make_set(setdir, [make_protein("P00001")])
            result = run_screening(setdir, Path(temp) / "screen")
            for name, ratio in result["summary"]["ratios"].items():
                self.assertIn("denominator", ratio, name)
                self.assertIn("denominator_definition", ratio, name)

    def test_figure_labelled_screening_not_success(self):
        with tempfile.TemporaryDirectory() as temp:
            setdir = Path(temp) / "set"
            make_set(setdir, [make_protein("P00001")])
            result = run_screening(setdir, Path(temp) / "screen", pilot=True)
            svg = (Path(result["run_dir"]) / "screening_overview.svg").read_text()
            self.assertIn("初筛建议", svg)
            self.assertIn("PILOT", svg)
            self.assertNotIn("成功率</text>", svg)
            data = read_json(Path(result["run_dir"]) / "screening_overview_data.json")
            self.assertEqual(data["not"], "实验成功率")

    def test_duplicate_and_scope_counting(self):
        dup_entry = {"entry_id": "entry-00003", "input": {"row": 3, "value": "P00001", "input_type": "accession"},
                     "identity_status": "resolved", "accession": "P00001", "gene": "SYNTHETIC_GENE",
                     "mapping": None, "fetch": None, "protein_file": "proteins/P00001.json",
                     "protein_sha256": "", "reference_selection": None, "processing_status": "ok",
                     "processing_error": "", "duplicate_of": "entry-00001",
                     "disposition": {"scope_status": "in_scope", "duplicate_input": True}}
        intra = make_protein("P00060", lambda p: p.update(topology="intracellular", location="unknown", features=[]))
        with tempfile.TemporaryDirectory() as temp:
            setdir = Path(temp) / "set"
            make_set(setdir, [make_protein("P00001"), intra], extra_entries=[dup_entry])
            # fix the duplicate entry hash after make_set wrote the protein file
            definition = read_json(setdir / "set_definition.json")
            definition["entries"][-1]["protein_sha256"] = file_hash(setdir / "proteins/P00001.json")
            write_json(setdir / "set_definition.json", definition)
            result = run_screening(setdir, Path(temp) / "screen")
            counts = result["summary"]["counts"]
            self.assertEqual(counts["input_rows"], 3)
            self.assertEqual(counts["duplicate_input_rows"], 1)
            self.assertEqual(counts["references_unique"], 2)  # duplicate not double-counted
            self.assertEqual(counts["out_of_scope"], 1)
            self.assertEqual(counts["evaluated_references"], 1)

    def test_legacy_assembly_branch_not_regressed(self):
        with tempfile.TemporaryDirectory() as temp:
            p = project()
            p["scaffold"] = scaffold()
            first = harness_run(p, temp)
            c = read_json(Path(first["run_dir"]) / "candidates.json")[0]
            p["reviews"] = [approve(c)]
            result = harness_run(p, temp)
            self.assertEqual(result["summary"]["fusion_sequences"], 1)


class CliChainTests(unittest.TestCase):
    def test_screen_cli_offline_chain(self):
        import ecd_snipr_cli as cli
        raw = Raw.entry("P00001", "CHAING")
        with tempfile.TemporaryDirectory() as temp:
            cache = Path(temp) / "cache"
            from ecd_snipr.uniprot import fetch
            self.assertEqual(fetch("P00001", cache, opener=lambda *a, **k: Response(raw))["status"], "fetched")
            target = Path(temp) / "targets.tsv"
            target.write_text("accession\nP00001\n", encoding="utf-8")
            code = cli.main(["screen", "--list", str(target), "--cache", str(cache),
                             "--outdir", str(Path(temp) / "out"), "--offline", "--pilot"])
            self.assertEqual(code, 0)
            runs = list((Path(temp) / "out").rglob("summary.json"))
            self.assertTrue(runs)
            summary = read_json(runs[0])
            self.assertEqual(summary["counts"]["evaluated_references"], 1)
            self.assertTrue(summary["pilot"])

    def test_build_set_cli_offline_partial(self):
        import ecd_snipr_cli as cli
        with tempfile.TemporaryDirectory() as temp:
            target = Path(temp) / "targets.tsv"
            target.write_text("accession\nP00009\n", encoding="utf-8")
            code = cli.main(["build-set", "--list", str(target), "--cache", str(Path(temp) / "cache"),
                             "--outdir", str(Path(temp) / "set"), "--offline"])
            self.assertEqual(code, 3)  # partial, not claimed complete
            definition = read_json(Path(temp) / "set/set_definition.json")
            self.assertEqual(definition["completeness"], "partial")
            self.assertEqual(definition["entries"][0]["identity_status"], "unresolved")


if __name__ == "__main__":
    unittest.main()
