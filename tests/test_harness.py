import copy
import json
import tempfile
import unittest
import urllib.error
import zipfile
from pathlib import Path
from unittest.mock import patch
from test_design import project, scaffold, approve, ROOT
from ecd_snipr.common import read_json
from ecd_snipr.contracts import EVIDENCE_FIELDS, evidence_v2, validate_project
from ecd_snipr.harness import run, verify_bundle
from ecd_snipr.ingest import inventory, map_rows
from ecd_snipr.uniprot import fetch, normalize
from test_assays import observation


class HarnessTests(unittest.TestCase):
    def test_end_to_end_cache_and_corruption_preservation(self):
        with tempfile.TemporaryDirectory() as temp:
            first = run(project(), temp)
            path = Path(first["run_dir"])
            self.assertTrue(verify_bundle(path))
            self.assertEqual(run(project(), temp, True)["execution"], "verified_cache_hit")
            (path / "candidate_plan.tsv").write_text("CORRUPTED")
            second = run(project(), temp, True)
            self.assertNotEqual(first["run_dir"], second["run_dir"])
            self.assertEqual((path / "candidate_plan.tsv").read_text(), "CORRUPTED")
            self.assertTrue(verify_bundle(second["run_dir"]))

    def test_no_overwrite_without_resume(self):
        with tempfile.TemporaryDirectory() as temp:
            run(project(), temp)
            with self.assertRaises(FileExistsError):
                run(project(), temp)

    def test_input_change_new_run(self):
        with tempfile.TemporaryDirectory() as temp:
            p = project()
            first = run(p, temp)
            p["rules"] = {"short_length": 40}
            self.assertNotEqual(first["run_dir"], run(p, temp)["run_dir"])

    def test_full_review_workflow(self):
        with tempfile.TemporaryDirectory() as temp:
            p = project()
            p["scaffold"] = scaffold()
            first = run(p, temp)
            c = read_json(Path(first["run_dir"]) / "candidates.json")[0]
            p["reviews"] = [approve(c)]
            result = run(p, temp)
            self.assertEqual(result["summary"]["fusion_sequences"], 1)
            self.assertTrue(verify_bundle(result["run_dir"]))

    def test_engine_change_invalidates_review(self):
        with tempfile.TemporaryDirectory() as temp:
            p = project()
            p["scaffold"] = scaffold()
            first = run(p, temp)
            p["reviews"] = [approve(read_json(Path(first["run_dir"]) / "candidates.json")[0])]
            with patch("ecd_snipr.harness.engine_hash", return_value="changed-engine"):
                self.assertEqual(run(p, temp)["summary"]["fusion_sequences"], 0)

    def test_construct_and_four_endpoint_pipeline(self):
        with tempfile.TemporaryDirectory() as temp:
            p = project()
            p["scaffold"] = scaffold()
            first = run(p, temp)
            c = read_json(Path(first["run_dir"]) / "candidates.json")[0]
            p["reviews"] = [approve(c)]
            second = run(p, temp)
            c = read_json(Path(second["run_dir"]) / "candidates.json")[0]
            p["existing_constructs"] = [{"construct_id": "C1", "protein_id": c["protein_id"], "start": c["start"], "end": c["end"], "antigen_sequence": c["sequence"], "fusion_sequence": c["fusion_sequence"], "scaffold_id": c["scaffold_id"], "scaffold_version": c["scaffold_version"], "source_ref": {"record": "synthetic://construct"}}]
            p["observations"] = []
            for i, endpoint in enumerate(("surface_expression", "recognition_retention", "basal_activation", "induced_response")):
                o = observation()
                o["observation_id"], o["endpoint"] = f"OBS{i}", endpoint
                if endpoint in {"induced_response", "recognition_retention"}:
                    o["context"].update(stimulus="binding_sender", sender_id="SYNTHETIC-SENDER", antibody_id="SYNTHETIC-ANTIBODY")
                p["observations"].append(o)
            result = run(p, temp)
            self.assertEqual(result["summary"]["assay_summary"]["eligible_measurements"], 4)
            c = read_json(Path(result["run_dir"]) / "candidates.json")[0]
            self.assertEqual(c["existing_construct_matches"][0]["scope"], "full_fusion_exact")
            self.assertEqual(c["functional_status"], "not_experimentally_validated")
            self.assertEqual(len(c["existing_construct_matches"][0]["observation_ids"]), 4)

    def test_duplicate_protein_ids_rejected(self):
        p = project()
        p["proteins"].append(copy.deepcopy(p["proteins"][0]))
        with self.assertRaises(ValueError):
            validate_project(p)

    def test_configuration_and_threshold_contract(self):
        for edits in ({"configuration": "antigen_sender"}, {"rules": {"short_length": 800, "long_length": 50}}):
            p = dict(project(), **edits)
            with self.assertRaises(ValueError):
                validate_project(p)

    def test_legacy_v2_preserved_and_extra_fields_rejected(self):
        r = dict.fromkeys(EVIDENCE_FIELDS, "")
        r.update(record_version="2.0", evidence_status="missing", source_authority_level="unknown", conflict_status="not_checked")
        evidence_v2(r)
        with tempfile.TemporaryDirectory() as temp:
            p = project()
            p["evidence_records_v2"] = [r]
            result = run(p, temp)
            self.assertEqual(json.loads((Path(result["run_dir"]) / "evidence_records_v2.jsonl").read_text()), r)
        r["construct_id"] = "cannot_extend_v2_silently"
        with self.assertRaises(ValueError):
            evidence_v2(r)

    def test_missing_legacy_evidence_not_zero(self):
        r = dict.fromkeys(EVIDENCE_FIELDS, "")
        r.update(record_version="2.0", evidence_status="missing", source_authority_level="unknown", conflict_status="not_checked", normalized_value="0")
        with self.assertRaises(ValueError):
            evidence_v2(r)


class ImportTests(unittest.TestCase):
    def test_tsv_mapping_preserves_y_n_blank_and_row(self):
        data = inventory(ROOT / "examples/synthetic.tsv")
        records = map_rows(data, read_json(ROOT / "examples/column_mapping.json"))
        self.assertEqual(records[0]["values"]["snipr_assay_raw"], "Y")
        self.assertEqual(records[0]["values"]["readout_raw"], "")
        self.assertEqual(records[1]["source"]["row"], 3)

    def test_wrong_sheet_is_error(self):
        data = inventory(ROOT / "examples/synthetic.tsv")
        mapping = read_json(ROOT / "examples/column_mapping.json")
        mapping["sheet"] = "synthetic "
        with self.assertRaises(ValueError):
            map_rows(data, mapping)

    def test_xlsx_preserves_formula_and_raw_cell(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "synthetic.xlsx"
            with zipfile.ZipFile(path, "w") as z:
                z.writestr("xl/workbook.xml", '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="test " sheetId="1" r:id="r1"/></sheets></workbook>')
                z.writestr("xl/_rels/workbook.xml.rels", '<Relationships><Relationship Id="r1" Target="worksheets/sheet1.xml"/></Relationships>')
                z.writestr("xl/worksheets/sheet1.xml", '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData><row r="2"><c r="A2" t="inlineStr"><is><t>Y</t></is></c><c r="C2"><f>1+1</f><v>2</v></c></row></sheetData><mergeCells><mergeCell ref="A2:B2"/></mergeCells></worksheet>')
            result = inventory(path)
            sheet = result["sheets"][0]
            self.assertEqual(sheet["name"], "test ")
            self.assertEqual(sheet["rows"][0]["cells"][1]["formula"], "1+1")
            self.assertEqual(sheet["merged_cells"], ["A2:B2"])


class UniProtTests(unittest.TestCase):
    def raw(self, description="Extracellular"):
        return {"primaryAccession": "P00001", "entryType": "UniProtKB reviewed (Swiss-Prot)", "sequence": {"value": "A" * 100}, "organism": {"taxonId": 9606}, "entryAudit": {"entryVersion": 1}, "features": [{"type": "Topological domain", "description": description, "location": {"start": {"value": 11}, "end": {"value": 70}}}, {"type": "Transmembrane", "location": {"start": {"value": 71}, "end": {"value": 90}}}]}

    def test_uniprot_preserves_isoform_uncertainty(self):
        p = normalize(self.raw())
        self.assertTrue(p["isoform_ambiguous"])
        self.assertEqual(p["topology"], "type_i")

    def test_lumenal_not_extracellular(self):
        p = normalize(self.raw("Lumenal"))
        self.assertEqual(p["features"][0]["kind"], "lumenal")
        self.assertEqual(p["topology"], "unknown")

    def test_eco_and_fuzzy_preserved(self):
        raw = self.raw()
        raw["features"][0]["location"]["start"]["modifier"] = "UNSURE"
        raw["features"][0]["evidences"] = [{"evidenceCode": "ECO:0000255"}]
        f = normalize(raw)["features"][0]
        self.assertEqual(f["boundary_status"], "fuzzy")
        self.assertEqual(f["evidence"]["kind"], "prediction")
        self.assertTrue(f["raw_feature"]["evidences"])

    def test_offline_no_network(self):
        with tempfile.TemporaryDirectory() as temp:
            with patch("urllib.request.urlopen", side_effect=AssertionError("network forbidden")):
                self.assertEqual(fetch("P00001", temp, offline=True)["status"], "missing")

    def test_retry_bounded(self):
        calls = []
        def fail(*args, **kwargs):
            calls.append(1)
            raise OSError("simulated network failure")
        with tempfile.TemporaryDirectory() as temp:
            self.assertEqual(fetch("P00001", temp, opener=fail, sleeper=lambda _: None)["status"], "error")
        self.assertEqual(len(calls), 3)

    def test_mock_fetch_and_cache_checksum(self):
        class Response:
            def __enter__(s): return s
            def __exit__(s, *args): return False
            def read(s): return json.dumps(self.raw()).encode()
        with tempfile.TemporaryDirectory() as temp:
            a = fetch("P00001", temp, opener=lambda *a, **k: Response())
            self.assertEqual(a["status"], "fetched")
            self.assertEqual(fetch("P00001", temp, offline=True)["status"], "cached")
            Path(a["cache_path"]).write_text("broken")
            self.assertEqual(fetch("P00001", temp, offline=True)["status"], "missing")

    def test_refresh_keeps_old_response_and_retrieval_record(self):
        raw = self.raw()
        class Response:
            def __enter__(s): return s
            def __exit__(s, *args): return False
            def read(s): return json.dumps(raw).encode()
        with tempfile.TemporaryDirectory() as temp:
            first = fetch("P00001", temp, opener=lambda *a, **k: Response())
            raw["entryAudit"]["entryVersion"] = 2
            second = fetch("P00001", temp, refresh=True, opener=lambda *a, **k: Response())
            self.assertNotEqual(first["cache_path"], second["cache_path"])
            self.assertTrue(Path(first["cache_path"]).is_file())
            self.assertEqual(len(list(Path(temp).glob("*.retrieval.json"))), 2)


if __name__ == "__main__":
    unittest.main()
