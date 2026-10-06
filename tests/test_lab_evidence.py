"""Synthetic regression only. No private workbook cells or target lists."""

import copy
import json
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from ecd_snipr.common import digest, file_hash, read_json, write_json
from ecd_snipr.contracts import evidence_v2
from ecd_snipr.harness import run, verify_bundle
from ecd_snipr.ingest import inventory, column_name
from ecd_snipr.lab_intake import extract, header_key
from ecd_snipr.lab_linkage import (assay_records, fragment_records, identity_conflicts, load_references,
                                   proposed_links, reviewed_assay_conflicts)
from ecd_snipr.lab_pipeline import run_lab_intake


def raw(rows):
    return {"source_name": "SYNTHETIC.xlsx", "source_sha256": "synthetic-hash", "sheets": [
        {"name": "Synthetic sheet ", "merged_cells": [], "rows": [
            {"row": i, "cells": [{"coordinate": f"{column_name(j)}{i}", "value": v, "formula": None}
                                  for j, v in enumerate(row, 1)]} for i, row in enumerate(rows, 1)]}]}


def fixture():
    return raw([["AgID", "gene", "Entry", "ECD_aa", "self-activation of Ag-SNIPR", "SNIPR Assay（Y/N）", "Ag-SNIPR plasmid", "note", "note"],
                ["SYN1", "SYN_GENE", "P12345", "ACDEFG", "0", "Y", "SYNPLASMID", "first", "second"],
                ["SYN2", "SYN_GENE", "P12345", "ACDEF", "", "N", "SYNPLASMID", "", ""]])


def semantics_for(rows):
    _, template = assay_records(rows)
    r = template[0]
    r.update(decision="confirmed", reviewer="synthetic_reviewer", reviewed_at="2026-01-01", rationale="synthetic definition")
    r["definition"] = {"endpoint": "basal_activation", "unit": "percent", "denominator": "synthetic live receiver cells",
                       "channel": "synthetic reporter", "gate_path": "live/receiver", "construct_id": "SYN_C1",
                       "context": {"batch_id": "SYN_B1", "host_cell": "synthetic", "sender_id": "none",
                                   "antibody_id": "not_applicable", "stimulus": "none", "condition_id": "SYN_COND",
                                   "timepoint": "SYN_T1", "replicate_id": "SYN_R1"}}
    return {"schema_version": "1.0", "source_sha256": "synthetic-hash", "reviews": [r]}


CONSTRUCTS = {"SYN_C1": {"audit_status": "sequence_consistent_not_functionally_validated", "antigen_sequence": "ACDEFG"}}


class LayoutTests(unittest.TestCase):
    def test_correct_and_legacy_typo_headers_unicode(self):
        a = extract(fixture())
        self.assertEqual(a["counts"]["nonblank_data_cells_by_field"]["self_activation"], 1)
        self.assertEqual(a["counts"]["nonblank_data_cells_by_field"]["receiver_plasmid"], 2)
        self.assertEqual(a["counts"]["nonblank_data_cells_by_field"]["snipr_assay_status"], 2)
        self.assertEqual(header_key(" SNIPR Assay（Y/N） "), header_key("SNIPR Assay(Y/N)"))
        f = fixture()
        f["sheets"][0]["rows"][0]["cells"][4]["value"] = "self-activation of Ag-SNPIR"
        self.assertEqual(extract(f)["counts"]["nonblank_data_cells_by_field"]["self_activation"], 1)

    def test_duplicate_headers_and_raw_values_preserved(self):
        e = extract(fixture())
        notes = [c["raw_value"] for c in e["rows"][0]["fields"] if c["field"] == "note"]
        self.assertEqual(notes, ["first", "second"])
        self.assertEqual(len(e["cells"]), 27)
        self.assertEqual(e["rows"][0]["sheet"], "Synthetic sheet ")

    def test_repeated_header_changed_columns(self):
        d = raw([["gene", "Entry", "ECD_aa"], ["SYN_A", "P12345", "ACDEFG"],
                 ["Entry", "gene", "ECD_nt"], ["P23456", "SYN_B", "GCTTGTGAT"]])
        e = extract(d)
        self.assertEqual(len(e["headers"]), 2)
        self.assertEqual(e["rows"][1]["fields"][0]["field"], "accession")
        self.assertEqual(len(e["rows"]), 2)

    def test_shifted_layout_quarantined(self):
        d = raw([["gene", "Entry", "ECD_aa"], ["SYN", "wrong", "P12345"]])
        rows = extract(d)["rows"]
        self.assertEqual(rows[0]["layout_status"], "suspect_or_composite_identity")
        self.assertFalse(fragment_records(rows, {})[0]["sequence"])

    def test_unknown_and_headerless_cells_retained(self):
        d = raw([["unknown", "private note"], ["other", "0"]])
        e = extract(d)
        self.assertEqual(len(e["cells"]), 4)
        self.assertEqual(e["rows"][0]["layout_status"], "no_header")

    def test_explicit_layout_hash_and_overlap(self):
        layout = {"source_sha256": "wrong", "regions": []}
        with self.assertRaises(ValueError):
            extract(fixture(), layout)
        region = {"sheet": "Synthetic sheet ", "start_row": 2, "end_row": 3,
                  "columns": {"A": "agid", "D": "ecd_aa"}, "rationale": "synthetic mapping"}
        layout = {"source_sha256": "synthetic-hash", "regions": [region, region]}
        with self.assertRaises(ValueError):
            extract(fixture(), layout)


class FragmentTests(unittest.TestCase):
    def setUp(self):
        self.rows = extract(fixture())["rows"]
        self.ref = {"P12345": {"protein_id": "SYN_REF", "accession": "P12345", "isoform": "P12345:canonical",
                               "sequence": "MACDEFGHIK", "topology": "type_i", "features": [], "evidence": {"source": "synthetic"}}}

    def test_exact_reference_not_lab_isoform_confirmation(self):
        f = fragment_records(self.rows, self.ref)[0]
        self.assertEqual(f["matches"][0]["start"], 2)
        self.assertEqual(f["matches"][0]["end"], 7)
        self.assertEqual(f["experimental_isoform_confirmation"], "not_performed")
        self.assertEqual(f["full_construct_identity"], "not_established")

    def test_ambiguous_positions_preserved(self):
        self.ref["P12345"]["sequence"] = "ACDEFGACDEFG"
        self.assertEqual(fragment_records(self.rows, self.ref)[0]["reference_status"], "ambiguous_multiple_positions")

    def test_core_invalid_sequence_does_not_stop_other_rows(self):
        self.rows[0]["fields"][3]["raw_value"] = "AXC*"
        f = fragment_records(self.rows, self.ref)
        self.assertFalse(f[0]["sequence"])
        self.assertTrue(f[1]["sequence"])

    def test_no_exact_match_not_wrong_construct_label(self):
        self.ref["P12345"]["sequence"] = "MMMMMMMMMM"
        self.assertEqual(fragment_records(self.rows, self.ref)[0]["reference_status"], "no_exact_match_not_proof_of_wrong_construct")

    def test_exclusions_cut_and_omitted_external_domains(self):
        self.ref["P12345"]["features"] = [
            {"kind": "signal_peptide", "start": 1, "end": 2},
            {"kind": "extracellular", "start": 2, "end": 9},
            {"kind": "domain", "start": 6, "end": 9},
            {"kind": "domain", "start": 8, "end": 9},
            {"kind": "domain", "start": 10, "end": 10}]
        m = fragment_records(self.rows, self.ref)[0]["matches"][0]
        self.assertEqual(len(m["retained_exclusions"]), 1)
        self.assertEqual(len(m["cut_domains"]), 1)
        self.assertEqual(len(m["omitted_extracellular_domains"]), 1)

    def test_nt_aa_conflict_preserved(self):
        rows = extract(raw([["gene", "Entry", "ECD_nt", "ECD_aa"], ["SYN", "P12345", "GCTTGTGAT", "ACDEFG"]]))["rows"]
        self.assertIn("ECD_nt_translation_vs_aa_mismatch", fragment_records(rows, self.ref)[0]["issues"])

    def test_untyped_ECD_kind_not_guessed_when_DNA_compatible(self):
        rows = extract(raw([["gene", "Entry", "ECD"], ["SYN", "P12345", "ACGT"]]))["rows"]
        self.assertFalse(fragment_records(rows, self.ref)[0]["sequence"])

    def test_shared_labels_variants_and_proposed_links_not_merge(self):
        f = fragment_records(self.rows, self.ref)
        self.assertEqual(len(identity_conflicts(f)), 1)
        edges, readiness = proposed_links(self.rows, f)
        self.assertTrue(all(e["status"] == "proposed_not_confirmed_construct_link" for e in edges))
        self.assertEqual(readiness[0]["candidate_unique_fragment_strings"], 2)

    def test_reference_checksum_error_isolated(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            write_json(root / "p.json", self.ref["P12345"])
            write_json(root / "set_definition.json", {"set_id": "synthetic", "completeness": "complete", "entries": [
                {"accession": "P12345", "protein_file": "p.json", "protein_sha256": "wrong"}]})
            refs, errors, _ = load_references(root, {"P12345"})
            self.assertFalse(refs)
            self.assertEqual(len(errors), 1)

    def test_plasmid_Y_and_checkmarks_are_not_shared_identity(self):
        for row in self.rows:
            row["fields"][6]["raw_value"] = "Y"
        f = fragment_records(self.rows, self.ref)
        self.assertFalse(identity_conflicts(f))
        edges, _ = proposed_links(self.rows, f)
        self.assertFalse(any("receiver_plasmid" in e["matching_fields"] for e in edges))


class SemanticTests(unittest.TestCase):
    def setUp(self):
        self.rows = extract(fixture())["rows"]

    def test_raw_zero_and_Y_N_and_blank_not_labels(self):
        records, _ = assay_records(self.rows)
        self.assertEqual(records[0]["raw_value"], "0")
        self.assertIsNone(records[0]["normalized_value"])
        self.assertEqual(records[0]["reported_endpoint"], "basal_activation")
        self.assertEqual(records[1]["record_kind"], "project_status_not_function")
        self.assertTrue(any(r["eligibility"] == "missing" for r in records))
        self.assertFalse(any(r["eligibility"].startswith("eligible") for r in records))

    def test_explicit_confirmed_definition_can_be_quantitative_not_success(self):
        o = assay_records(self.rows, semantics_for(self.rows), CONSTRUCTS)[0][0]
        self.assertEqual(o["normalized_value"], 0)
        self.assertEqual(o["eligibility"], "eligible_measurement_not_success_label")

    def test_missing_actual_construct_still_unresolved(self):
        o = assay_records(self.rows, semantics_for(self.rows))[0][0]
        self.assertIsNone(o["normalized_value"])
        self.assertIn("construct_identity_unresolved", o["issues"])

    def test_explicit_review_cannot_transfer_to_different_fragment(self):
        c = copy.deepcopy(CONSTRUCTS)
        c["SYN_C1"]["antigen_sequence"] = "AAAAA"
        o = assay_records(self.rows, semantics_for(self.rows), c)[0][0]
        self.assertIsNone(o["normalized_value"])
        self.assertIn("reviewed_construct_conflicts_with_row_fragment", o["issues"])

    def test_stale_hash_and_unknown_cell_rejected(self):
        for key, value in (("raw_record_digest", "stale"), ("cell_id", "absent")):
            s = semantics_for(self.rows)
            s["reviews"][0][key] = value
            with self.assertRaises(ValueError):
                assay_records(self.rows, s, CONSTRUCTS)

    def test_review_cannot_override_raw_value(self):
        s = semantics_for(self.rows)
        s["reviews"][0]["definition"]["raw_value"] = "99"
        with self.assertRaises(ValueError):
            assay_records(self.rows, s, CONSTRUCTS)

    def test_Y_N_cannot_be_promoted_by_review(self):
        s = semantics_for(self.rows)
        c = next(c for c in self.rows[0]["fields"] if c["field"] == "snipr_assay_status")
        s["reviews"][0].update(cell_id=c["cell_id"], raw_record_digest=digest(c))
        with self.assertRaises(ValueError):
            assay_records(self.rows, s, CONSTRUCTS)

    def test_formula_requires_cache_review(self):
        self.rows[0]["fields"][4]["raw_formula"] = "1-1"
        s = semantics_for(self.rows)
        with self.assertRaises(ValueError):
            assay_records(self.rows, s, CONSTRUCTS)
        s["reviews"][0]["formula_cache_reviewed"] = True
        self.assertEqual(assay_records(self.rows, s, CONSTRUCTS)[0][0]["normalized_value"], 0)

    def test_shared_formula_follower_requires_review_even_without_text(self):
        cell = self.rows[0]["fields"][4]
        cell.update(raw_formula=None, formula_present=True, formula_attributes={"t": "shared", "si": "0"},
                    shared_formula_anchors=[{"coordinate": "E1", "formula": "A1*0"}])
        with self.assertRaises(ValueError):
            assay_records(self.rows, semantics_for(self.rows), CONSTRUCTS)

    def test_cached_fragment_needs_sequence_cache_review(self):
        self.rows[0]["fields"][3]["formula_present"] = True
        s = semantics_for(self.rows)
        o = assay_records(self.rows, s, CONSTRUCTS)[0][0]
        self.assertIn("sequence_formula_cache_requires_review", o["issues"])
        s["reviews"][0]["sequence_formula_cache_reviewed"] = True
        self.assertEqual(assay_records(self.rows, s, CONSTRUCTS)[0][0]["eligibility"], "eligible_measurement_not_success_label")

    def test_unresolved_contexts_not_spurious_conflicts(self):
        o = assay_records(self.rows)[0]
        self.assertEqual(reviewed_assay_conflicts(o), [])

    def test_batches_not_pooled(self):
        o = assay_records(self.rows, semantics_for(self.rows), CONSTRUCTS)[0][0]
        b = copy.deepcopy(o)
        b["observation_id"], b["raw_value"], b["context"]["batch_id"] = "synthetic2", "20", "otherbatch"
        self.assertEqual(reviewed_assay_conflicts([o, b]), [])
        b["context"]["batch_id"] = o["context"]["batch_id"]
        self.assertEqual(len(reviewed_assay_conflicts([o, b])), 1)


class PipelineTests(unittest.TestCase):
    def test_review_and_real_audit_branch_end_to_end(self):
        project = read_json(Path(__file__).resolve().parents[1] / "examples/synthetic_project.json")
        p = project["proteins"][0]
        antigen = p["sequence"][10:20]
        project["existing_constructs"] = [{"construct_id": "SYN_C1", "protein_id": p["protein_id"],
            "start": 11, "end": 20, "antigen_sequence": antigen, "fusion_sequence": "M" + antigen + "GGS",
            "scaffold_id": "SYN_SCAFFOLD", "scaffold_version": "1", "source_ref": {"record": "synthetic://record"}}]
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            construct_run = run(project, root / "construct")
            source = root / "synthetic.tsv"
            source.write_text("gene\tEntry\tECD_aa\tself-activation of Ag-SNIPR\nSYN\tP12345\t" + antigen + "\t0\n")
            rows = extract(inventory(source))["rows"]
            sem = semantics_for(rows)
            sem["source_sha256"] = file_hash(source)
            write_json(root / "review.json", sem)
            result = run_lab_intake(source, root / "lab", semantics_path=root / "review.json", construct_run=construct_run["run_dir"])
            self.assertEqual(result["summary"]["counts"]["eligible_measurements"], 1)
            self.assertTrue(verify_bundle(result["run_dir"]))

    def test_generated_lab_artifact_names_excluded_from_package(self):
        from manage_skill import forbidden
        for name in ("examples/raw_inventory.json", "tests/source_cells.tsv", "lab_evidence_output/result.json"):
            self.assertTrue(forbidden(Path(name)))

    def test_end_to_end_no_scaffold_source_unchanged_resume_and_tamper(self):
        source = Path(__file__).resolve().parents[1] / "examples/synthetic_lab.tsv"
        before = file_hash(source)
        with tempfile.TemporaryDirectory() as temp:
            result = run_lab_intake(source, temp)
            root = Path(result["run_dir"])
            self.assertEqual(result["summary"]["counts"]["explicit_self_activation_source_cells"], 1)
            self.assertTrue(verify_bundle(root))
            self.assertEqual(run_lab_intake(source, temp, resume=True)["execution"], "reused_verified")
            for line in (root / "evidence_records_v2.jsonl").read_text().splitlines():
                evidence_v2(json.loads(line))
            (root / "summary.json").write_text("{}")
            with self.assertRaises(FileExistsError):
                run_lab_intake(source, temp, resume=True)
        self.assertEqual(before, file_hash(source))

    def test_xlsx_original_formula_merge_unicode_and_duplicate_headers(self):
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / "synthetic.xlsx"
            ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
            with zipfile.ZipFile(source, "w") as z:
                z.writestr("xl/workbook.xml", f'<workbook xmlns="{ns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Synthetic sheet " sheetId="1" r:id="rId1"/></sheets></workbook>')
                z.writestr("xl/_rels/workbook.xml.rels", '<Relationships><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
                root = ET.Element("worksheet", xmlns=ns)
                data = ET.SubElement(root, "sheetData")
                for row in fixture()["sheets"][0]["rows"]:
                    node = ET.SubElement(data, "row", r=str(row["row"]))
                    for cell in row["cells"]:
                        n = ET.SubElement(node, "c", r=cell["coordinate"], t="inlineStr")
                        ET.SubElement(ET.SubElement(n, "is"), "t").text = cell["value"]
                        if cell["coordinate"] == "E2":
                            ET.SubElement(n, "f", t="shared", si="0", ref="E2:E3").text = "1-1"
                        if cell["coordinate"] == "E3":
                            ET.SubElement(n, "f", t="shared", si="0")
                ET.SubElement(ET.SubElement(root, "mergeCells"), "mergeCell", ref="H3:I3")
                z.writestr("xl/worksheets/sheet1.xml", ET.tostring(root))
            before = file_hash(source)
            result = run_lab_intake(source, Path(temp) / "output")
            inv = read_json(Path(result["run_dir"]) / "raw_inventory.json")
            self.assertEqual(inv["sheets"][0]["merged_cells"], ["H3:I3"])
            self.assertEqual(inv["sheets"][0]["rows"][1]["cells"][4]["formula"], "1-1")
            follower = inv["sheets"][0]["rows"][2]["cells"][4]
            self.assertTrue(follower["formula_present"])
            self.assertIsNone(follower["formula"])
            self.assertEqual(follower["shared_formula_anchors"][0]["coordinate"], "E2")
            self.assertEqual(follower["formula_attributes"]["si"], "0")
            self.assertEqual(before, file_hash(source))


if __name__ == "__main__":
    unittest.main()
