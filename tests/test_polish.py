"""Regression tests for the v0.3.1 polish round.

Covers: gene-mapping pagination (no false 'resolved' from a truncated first
page), per-page retry with bounded Retry-After on search/fetch, per-entry
isolation of corrupted cached protein records (never a whole-run abort),
resolved-only gene/reference counts in summary.json, and isoform inventory /
isoform-difference recording that never changes a screening class.
"""
import io
import json
import tempfile
import unittest
import urllib.error
from pathlib import Path
from test_design import ROOT
from test_screening import Raw, Response, make_protein, make_set
from ecd_snipr.acquisition import map_gene_name, search_accessions, select_analysis_reference
from ecd_snipr.common import read_json
from ecd_snipr.harness import verify_bundle
from ecd_snipr.screening import run_screening
from ecd_snipr.uniprot import fetch, normalize, retry_delay


def http_error(code, headers=None):
    return urllib.error.HTTPError("https://rest.uniprot.org/x", code, "err", headers or {}, io.BytesIO(b""))


class GeneMappingPaginationTests(unittest.TestCase):
    def test_gene_mapping_follows_pagination_before_resolving(self):
        # Page 1 yields one exact-gene hit; page 2 yields a second. A truncated
        # single-page read would falsely report "resolved".
        def opener(request, **kwargs):
            url = request.full_url
            if "cursor=PAGE2" in url:
                return Response({"results": [Raw.entry("P22222", "SPLITG")]}, {"X-UniProt-Release": "2026_02"})
            return Response({"results": [Raw.entry("P11111", "SPLITG")]},
                            {"Link": '<https://rest.uniprot.org/uniprotkb/search?query=x&cursor=PAGE2>; rel="next"',
                             "X-UniProt-Release": "2026_02"})
        with tempfile.TemporaryDirectory() as temp:
            mapping = map_gene_name("SPLITG", temp, opener=opener, sleeper=lambda _: None)
            self.assertEqual(mapping["status"], "ambiguous")
            self.assertEqual({c["accession"] for c in mapping["candidates"]}, {"P11111", "P22222"})
            self.assertEqual(mapping["queries"][0]["uniprot_release"], "2026_02")

    def test_gene_mapping_pagination_cap_is_error_not_silent_truncation(self):
        def opener(request, **kwargs):
            return Response({"results": [Raw.entry("P33333", "HUGE")]},
                            {"Link": '<https://rest.uniprot.org/uniprotkb/search?query=x&cursor=NEXT>; rel="next"'})
        with tempfile.TemporaryDirectory() as temp:
            from ecd_snipr.acquisition import _query_with_retry
            hits, release, error = _query_with_retry("gene_exact:HUGE", opener, lambda _: None, max_pages=2)
            self.assertEqual(error, "pagination_cap_reached_incomplete_candidate_list")
            mapping = map_gene_name("HUGE", temp, opener=opener, sleeper=lambda _: None)
            self.assertEqual(mapping["status"], "error")  # honest failure, not a guessed resolution

    def test_search_accessions_retries_transient_page_error(self):
        calls = {"n": 0}
        def opener(request, **kwargs):
            calls["n"] += 1
            if calls["n"] == 2:
                raise http_error(503)  # transient failure on page 2, retried
            page = {"results": [Raw.entry(f"P{calls['n']:05d}", "G")]}
            headers = {"Link": '<https://rest.uniprot.org/x?query=q&cursor=ABC>; rel="next"'} if calls["n"] == 1 else {}
            return Response(page, headers)
        result = search_accessions("(organism_id:9606)", opener=opener, sleeper=lambda _: None)
        self.assertEqual(result["status"], "fetched")
        self.assertEqual(len(result["entries"]), 2)
        self.assertEqual(result["completeness"], "complete")


class RetryAfterTests(unittest.TestCase):
    def test_retry_delay_honors_bounded_retry_after(self):
        self.assertEqual(retry_delay(http_error(429, {"Retry-After": "7"}), 0), 7)
        self.assertEqual(retry_delay(http_error(503, {"Retry-After": "3600"}), 0), 60)  # capped
        self.assertEqual(retry_delay(http_error(429), 1), 2)  # no header: exponential fallback
        self.assertEqual(retry_delay(OSError("network"), 0), 1)

    def test_fetch_429_retry_after_and_404_no_retry(self):
        slept = []
        calls = {"n": 0}
        def opener(request, **kwargs):
            calls["n"] += 1
            if calls["n"] == 1:
                raise http_error(429, {"Retry-After": "5"})
            return Response(Raw.entry("P00001"))
        with tempfile.TemporaryDirectory() as temp:
            status = fetch("P00001", temp, opener=opener, sleeper=slept.append)
            self.assertEqual(status["status"], "fetched")
            self.assertEqual(slept, [5])
            def opener_404(request, **kwargs):
                raise http_error(404)
            status = fetch("P00009", temp, opener=opener_404, sleeper=slept.append)
            self.assertEqual(status["status"], "error")
            self.assertEqual(slept, [5])  # non-retryable: no additional sleep


class CorruptRecordIsolationTests(unittest.TestCase):
    def test_corrupt_protein_record_isolated_and_marks_partial(self):
        with tempfile.TemporaryDirectory() as temp:
            setdir = Path(temp) / "set"
            make_set(setdir, [make_protein("P00001"), make_protein("P00002")])
            bad = setdir / "proteins" / "P00002.json"
            bad.write_text(bad.read_text() + "tampered")
            result = run_screening(setdir, Path(temp) / "screen")
            records = read_json(Path(result["run_dir"]) / "screening.json")
            ok = next(r for r in records if r["accession"] == "P00001")
            broken = next(r for r in records if r["accession"] == "P00002")
            self.assertEqual(ok["screening_recommendation"], "standard_candidate")
            self.assertEqual(broken["processing_status"], "technical_failure")
            self.assertIsNone(broken["screening_recommendation"])
            self.assertIn("checksum mismatch", broken["processing_error"])
            # The acquisition definition said complete; the screen-stage failure
            # must still mark the run partial.
            self.assertEqual(result["summary"]["completeness"]["state"], "partial")
            self.assertEqual(result["summary"]["counts"]["technical_failures"], 1)
            self.assertTrue(verify_bundle(result["run_dir"]))


class SummaryCountingTests(unittest.TestCase):
    def test_summary_counts_only_resolved_genes_and_references(self):
        failed_entry = {"entry_id": "entry-00002", "input": {"row": 2, "value": "P99999", "input_type": "accession"},
                        "identity_status": "unresolved", "accession": "P99999", "gene": "GHOST",
                        "mapping": None, "fetch": {"status": "error", "reason": "simulated"},
                        "protein_file": "", "protein_sha256": "", "reference_selection": None,
                        "processing_status": "failed", "processing_error": "OSError: simulated",
                        "disposition": {"scope_status": "unevaluated", "note": "technical failure"}}
        with tempfile.TemporaryDirectory() as temp:
            setdir = Path(temp) / "set"
            make_set(setdir, [make_protein("P00001")], extra_entries=[failed_entry], failed_rows=[failed_entry])
            result = run_screening(setdir, Path(temp) / "screen")
            counts = result["summary"]["counts"]
            self.assertEqual(counts["references_unique"], 1)  # not 2: unresolved row is not a reference
            self.assertEqual(counts["genes_unique"], 1)       # not 2: GHOST never resolved
            self.assertEqual(counts["unresolved_identities"], 1)
            self.assertEqual(counts["technical_failures"], 1)


def lag3_like_raw():
    """Reviewed entry with two annotated isoforms and an in-ECD difference."""
    raw = Raw.entry("P18627", "LAG3X")
    raw["comments"].append({"commentType": "ALTERNATIVE PRODUCTS", "events": ["Alternative splicing"],
                            "isoforms": [{"name": {"value": "1"}, "isoformIds": ["P18627-1"], "isoformSequenceStatus": "Displayed"},
                                         {"name": {"value": "2"}, "isoformIds": ["P18627-2"], "isoformSequenceStatus": "Described"}]})
    raw["features"].append({"type": "Alternative sequence", "description": "in isoform 2",
                            "location": {"start": {"value": 20, "modifier": "EXACT"},
                                         "end": {"value": 30, "modifier": "EXACT"}},
                            "featureId": "VSP_000001",
                            "alternativeSequence": {"originalSequence": "VTPKSFGS", "alternativeSequences": ["GQPQVGKE"]}})
    return raw


class IsoformRecordingTests(unittest.TestCase):
    def test_isoform_inventory_recorded_in_normalize_and_reference_selection(self):
        p = select_analysis_reference(normalize(lag3_like_raw()))
        ap = p["alternative_products"]
        self.assertEqual(ap["isoform_count"], 2)
        self.assertEqual([i["name"] for i in ap["isoforms"]], ["1", "2"])
        self.assertEqual(ap["events"], ["Alternative splicing"])
        self.assertEqual(p["isoform_differences"][0]["feature_id"], "VSP_000001")
        sel = p["reference_selection"]
        self.assertEqual(sel["annotated_isoform_count"], 2)
        self.assertEqual(sel["isoform_comparison"], "not_evaluated")
        self.assertIn("2 个 isoform", sel["isoform_comparison_note"])
        # No ALTERNATIVE PRODUCTS comment: recorded as not annotated, not as zero.
        plain = select_analysis_reference(normalize(Raw.entry("P00001")))
        self.assertIsNone(plain["alternative_products"])
        self.assertIsNone(plain["reference_selection"]["annotated_isoform_count"])
        self.assertIn("不证明不存在", plain["reference_selection"]["isoform_comparison_note"])

    def test_isoform_difference_within_candidate_recorded_never_downgrades(self):
        with tempfile.TemporaryDirectory() as temp:
            setdir = Path(temp) / "set"
            make_set(setdir, [select_analysis_reference(normalize(lag3_like_raw()))])
            result = run_screening(setdir, Path(temp) / "screen")
            rec = read_json(Path(result["run_dir"]) / "screening.json")[0]
            self.assertEqual(rec["screening_recommendation"], "standard_candidate")  # LAG3 deferral preserved
            diffs = rec["evidence"]["isoform_differences_within_candidate"]
            self.assertEqual(diffs[0]["feature_id"], "VSP_000001")
            self.assertTrue(any("可变序列差异落在候选区间内" in m for m in rec["missing_info"]))
            self.assertTrue(any("2 个 isoform" in m for m in rec["missing_info"]))

    def test_isoform_difference_outside_candidate_not_surfaced(self):
        raw = lag3_like_raw()
        raw["features"][-1]["location"] = {"start": {"value": 95, "modifier": "EXACT"},
                                           "end": {"value": 98, "modifier": "EXACT"}}  # cytoplasmic tail
        with tempfile.TemporaryDirectory() as temp:
            setdir = Path(temp) / "set"
            make_set(setdir, [select_analysis_reference(normalize(raw))])
            result = run_screening(setdir, Path(temp) / "screen")
            rec = read_json(Path(result["run_dir"]) / "screening.json")[0]
            self.assertNotIn("isoform_differences_within_candidate", rec["evidence"])
            self.assertFalse(any("可变序列差异落在候选区间内" in m for m in rec["missing_info"]))


if __name__ == "__main__":
    unittest.main()
