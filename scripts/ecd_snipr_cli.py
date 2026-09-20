#!/usr/bin/env python3
"""Public CLI. Files are never overwritten by one-off import/export commands."""

import argparse
import json
import sys
from pathlib import Path
from ecd_snipr.acquisition import DEFAULT_SET_QUERY, build_set, read_target_list
from ecd_snipr.common import now, read_json, write_json
from ecd_snipr.harness import run, verify_bundle
from ecd_snipr.ingest import inventory, map_rows
from ecd_snipr.screening import run_screening
from ecd_snipr.uniprot import fetch, normalize


def new_output(path, value):
    if Path(path).exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")
    write_json(path, value)


def parser():
    p = argparse.ArgumentParser(description="Traceable antigen-ECD -> lab SNIPR receiver candidate design. No success prediction.")
    commands = p.add_subparsers(dest="command", required=True)
    for name in ("run", "smoke"):
        sub = commands.add_parser(name, help="Offline synthetic fixture pipeline" if name == "smoke" else "Execute a normalized project")
        if name == "run":
            sub.add_argument("--project", required=True)
            sub.add_argument("--reviews", help="Human-authored scoped review JSON list")
        sub.add_argument("--outdir", required=True)
        sub.add_argument("--resume", action="store_true")
    sub = commands.add_parser("inspect", help="Read XLSX/CSV/TSV cells without interpreting assay semantics")
    sub.add_argument("--input", required=True)
    sub.add_argument("--output", required=True)
    sub = commands.add_parser("map", help="Explicit column mapping; retains raw cell provenance")
    sub.add_argument("--inventory", required=True)
    sub.add_argument("--mapping", required=True)
    sub.add_argument("--output", required=True)
    sub = commands.add_parser("normalize-uniprot", help="Import an existing UniProt JSON; does not choose the experimental isoform")
    sub.add_argument("--input", required=True)
    sub.add_argument("--output", required=True)
    sub = commands.add_parser("fetch-uniprot", help="Optional accession-only network operation; never part of offline CI")
    sub.add_argument("--accession", required=True)
    sub.add_argument("--cache", required=True)
    sub.add_argument("--offline", action="store_true")
    sub.add_argument("--refresh", action="store_true")
    sub.add_argument("--output", required=True, help="Status record, including missing/error")
    sub = commands.add_parser("review-template", help="Generate pending reviews, never auto-approve")
    sub.add_argument("--run-dir", required=True)
    sub.add_argument("--output", required=True)
    sub = commands.add_parser("verify-run", help="Check every artifact checksum")
    sub.add_argument("--run-dir", required=True)
    sub = commands.add_parser("build-set", help="Build a membrane-protein analysis set from a target list or a public UniProt query")
    source = sub.add_mutually_exclusive_group(required=True)
    source.add_argument("--list", help="User target list (TSV/CSV/TXT): accessions or gene names; every row preserved")
    source.add_argument("--query", help="Public UniProt query, e.g. the documented human membrane default")
    sub.add_argument("--cache", required=True, help="Fetch/mapping cache directory (content-addressed, resumable)")
    sub.add_argument("--outdir", required=True, help="Set directory for set_definition.json and normalized proteins")
    sub.add_argument("--offline", action="store_true", help="Use cache only; no network")
    sub.add_argument("--refresh", action="store_true", help="Re-fetch even when a valid cache entry exists")
    sub.add_argument("--limit", type=int, help="Maximum query hits (marks the set truncated/partial)")
    sub = commands.add_parser("screen", help="Batch screening: set -> normalize -> recommendations -> summary. No per-protein JSON needed")
    source = sub.add_mutually_exclusive_group(required=True)
    source.add_argument("--set", help="Set directory containing set_definition.json from build-set")
    source.add_argument("--list", help="Target list; the set is built first, then screened (needs --cache)")
    source.add_argument("--query", help="Public UniProt query; the set is built first, then screened (needs --cache)")
    sub.add_argument("--cache", help="Fetch/mapping cache for the chained build-set step")
    sub.add_argument("--outdir", required=True, help="Screening output root (immutable run bundle per input)")
    sub.add_argument("--set-outdir", help="Where to write the chained set (default: <outdir>/analysis_set)")
    sub.add_argument("--lab-rules", help="Optional sourced lab experience rules JSON; absent rules are reported as 尚未纳入")
    sub.add_argument("--domain-policy", choices=["auto", "full_ecd_only", "full_ecd_and_domains"], default="auto",
                     help="auto: multi-pass references explicitly opt into sourced external-domain alternates")
    sub.add_argument("--pilot", action="store_true", help="Label all outputs and the figure as a small pilot run")
    sub.add_argument("--review-per-stratum", type=int, default=2,
                     help="Pending manual-review rows per topology/class/evidence stratum (1-100, default 2); not an accuracy sample")
    sub.add_argument("--offline", action="store_true")
    sub.add_argument("--refresh", action="store_true")
    sub.add_argument("--limit", type=int)
    sub.add_argument("--resume", action="store_true")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    if args.command in {"run", "smoke"}:
        root = Path(__file__).resolve().parents[1]
        project = read_json(root / "examples/synthetic_project.json" if args.command == "smoke" else args.project)
        if getattr(args, "reviews", None):
            project["reviews"] = read_json(args.reviews)
        result = run(project, args.outdir, args.resume)
        if args.command == "smoke":
            summary = result["summary"]
            if summary["proteins"] != 1 or summary["candidates"] != 3 or summary["fusion_sequences"] != 0:
                raise ValueError("Fixture smoke invariant failed")
            if not verify_bundle(result["run_dir"]):
                raise ValueError("Fixture artifact integrity failed")
        print(json.dumps(result, ensure_ascii=False, indent=2))
    elif args.command == "inspect":
        data = inventory(args.input)
        new_output(args.output, data)
        print(json.dumps({"source_sha256": data["source_sha256"], "sheets": [{"name": s["name"], "rows": len(s["rows"])} for s in data["sheets"]]}, ensure_ascii=False))
    elif args.command == "map":
        result = map_rows(read_json(args.inventory), read_json(args.mapping))
        new_output(args.output, result)
        print(json.dumps({"mapped_rows": len(result), "status": "uninterpreted"}))
    elif args.command == "normalize-uniprot":
        new_output(args.output, normalize(read_json(args.input)))
        print("Imported; confirm actual isoform and source coordinates before design")
    elif args.command == "fetch-uniprot":
        if Path(args.output).exists():
            raise FileExistsError(args.output)
        result = fetch(args.accession, args.cache, args.offline, args.refresh)
        new_output(args.output, result)
        print(json.dumps(result))
        return 0 if result["status"] in {"fetched", "cached"} else 2
    elif args.command == "review-template":
        if not verify_bundle(args.run_dir):
            raise ValueError("Run integrity check failed")
        reviews = [{"review_key": c["review_key"], "candidate_id": c["candidate_id"], "decision": "pending", "reviewer": "", "reviewed_at": "", "rationale": "", "acknowledged_risks": [], "required_risks": sorted({f["code"] for f in c["risks"] if f["severity"] == "review"})} for c in read_json(Path(args.run_dir) / "candidates.json") if c["review_key"]]
        new_output(args.output, reviews)
        print("Pending review template only; authorized reviewer must fill decision and risk acknowledgements")
    elif args.command == "verify-run":
        verified = verify_bundle(args.run_dir)
        print(json.dumps({"verified": verified, "checked_at": now()}))
        return 0 if verified else 1
    elif args.command == "build-set":
        definition = _build_set_command(args)
        print(json.dumps({"set_id": definition["set_id"], "counts": definition["counts"],
                          "completeness": definition["completeness"]}, ensure_ascii=False, indent=2))
        return 0 if definition["completeness"] == "complete" else 3
    elif args.command == "screen":
        set_dir = args.set
        if not set_dir:
            if not args.cache:
                raise ValueError("--list/--query screening requires --cache for the chained build-set step")
            set_dir = args.set_outdir or str(Path(args.outdir) / "analysis_set")
            _build_set_command(args, set_dir)
        result = run_screening(set_dir, args.outdir, lab_rules_path=args.lab_rules,
                               domain_policy=args.domain_policy, pilot=args.pilot, resume=args.resume,
                               review_per_stratum=args.review_per_stratum)
        summary = result["summary"]
        print(json.dumps({"run_dir": result["run_dir"], "execution": result["execution"],
                          "set_completeness": summary["completeness"]["state"],
                          "counts": summary["counts"]}, ensure_ascii=False, indent=2))
        return 0 if summary["completeness"]["state"] == "complete" else 3
    return 0


def _build_set_command(args, outdir=None):
    rows = read_target_list(args.list) if getattr(args, "list", None) else None
    query = getattr(args, "query", None)
    if query == "default":
        query = DEFAULT_SET_QUERY
    return build_set(rows=rows, query=query, cache=args.cache,
                     outdir=outdir or args.outdir, offline=args.offline,
                     refresh=getattr(args, "refresh", False), limit=getattr(args, "limit", None))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, KeyError, TypeError, OSError) as exc:
        print(json.dumps({"status": "error", "reason": str(exc)}, ensure_ascii=False), file=sys.stderr)
        raise SystemExit(1)
