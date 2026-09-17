"""Immutable run bundles, content-addressed cache, state and integrity checks."""

from collections import Counter
from pathlib import Path
from . import __version__
from .common import canonical, digest, file_hash, now, read_json, tsv, write_json
from .contracts import validate_artifacts, validate_project
from .design import audit_construct, propose
from .observations import assay_summary, normalize_observation
from .provenance import claims_for

CANDIDATE_COLUMNS = ["candidate_id", "protein_id", "gene", "accession", "isoform", "topology", "antigen_form_type", "start", "end", "length", "design_status", "assembly_status", "risk_review_status", "rationale", "domains_retained", "domains_cut", "domains_omitted", "epitope_review", "risks", "scaffold_issues", "existing_construct_matches", "review_key", "synthetic_only"]


def engine_hash():
    root = Path(__file__).resolve().parents[1]
    return digest({str(p.relative_to(root)): file_hash(p) for p in sorted(root.rglob("*.py"))})


def verify_bundle(path):
    path = Path(path)
    try:
        manifest = read_json(path / "manifest.json")
        if manifest["state"] != "complete":
            return False
        return bool(manifest["artifacts"]) and all(Path(name).name == name and (path / name).is_file() and file_hash(path / name) == checksum for name, checksum in manifest["artifacts"].items())
    except (OSError, KeyError, ValueError):
        return False


def run(project, outdir, resume=False):
    validate_project(project)
    code_hash = engine_hash()
    run_key = digest({"project": project, "engine_sha256": code_hash})
    root = Path(outdir)
    root.mkdir(parents=True, exist_ok=True)
    base = root / run_key[:20]
    if resume and verify_bundle(base):
        return {"run_dir": str(base.resolve()), "execution": "verified_cache_hit", "summary": read_json(base / "summary.json")}
    run_dir = base
    attempt = 1
    while run_dir.exists():
        if not resume:
            raise FileExistsError("Run already exists; use --resume to verify cached artifacts, or a new output root")
        attempt += 1
        run_dir = root / f"{run_key[:20]}-attempt-{attempt}"
        if verify_bundle(run_dir):
            return {"run_dir": str(run_dir.resolve()), "execution": "verified_cache_hit", "summary": read_json(run_dir / "summary.json")}
    run_dir.mkdir()
    state = {"run_key": run_key, "engine_sha256": code_hash, "version": __version__, "started_at": now(), "state": "running", "phases": []}

    def event(phase, status, **details):
        item = dict(phase=phase, status=status, timestamp=now(), **details)
        with (run_dir / "events.jsonl").open("a", encoding="utf-8") as handle:
            handle.write(canonical(item) + "\n")
        state["phases"].append(item)
        write_json(run_dir / "state.json", state)

    try:
        write_json(run_dir / "input_snapshot.json", project)
        event("input_validation", "complete")
        candidates, diagnostics = [], []
        rules = dict(project.get("rules", {}), engine_sha256=code_hash)
        for protein in project["proteins"]:
            result, issues = propose(protein, project.get("scaffold"), rules, project.get("reviews", []))
            candidates.extend(result)
            diagnostics.append({"protein_id": protein["protein_id"], "candidate_count": len(result), "diagnostics": issues})
        event("candidate_design", "complete", proteins=len(diagnostics), candidates=len(candidates))
        proteins = {p["protein_id"]: p for p in project["proteins"]}
        constructs = {c["construct_id"]: audit_construct(c, proteins) for c in project.get("existing_constructs", [])}
        observations = [normalize_observation(o, constructs) for o in project.get("observations", [])]
        for c in candidates:
            c["existing_construct_matches"] = []
            for old in constructs.values():
                if (old.get("protein_id"), old.get("start"), old.get("end"), old.get("antigen_sequence")) == (c["protein_id"], c.get("start"), c.get("end"), c["sequence"]):
                    exact = bool(c["fusion_sequence"] and c["fusion_sequence"] == old.get("fusion_sequence") and c.get("scaffold_id") == old.get("scaffold_id") and c.get("scaffold_version") == old.get("scaffold_version"))
                    c["existing_construct_matches"].append({"construct_id": old["construct_id"], "scope": "full_fusion_exact" if exact else "antigen_only_not_equivalent_receptor", "audit_status": old["audit_status"], "observation_ids": [o["observation_id"] for o in observations if o.get("construct_id") == old["construct_id"]]})
        event("experiment_linkage", "complete", observations=len(observations))
        claims, links, conflicts = claims_for(project, candidates, diagnostics, observations)
        validate_artifacts(candidates, observations)
        event("artifact_validation", "complete")
        summary = {"configuration": project["configuration"], "proteins": len(proteins), "candidates": len(candidates),
                   "proteins_with_candidates": sum(d["candidate_count"] > 0 for d in diagnostics),
                   "design_status": dict(Counter(c["design_status"] for c in candidates)),
                   "fusion_sequences": sum(bool(c["fusion_sequence"]) for c in candidates),
                   "functional_success_probability": "not_computed", "conflicts": len(conflicts),
                   "assay_summary": assay_summary(observations, constructs),
                   "interpretation": "Computational workflow completed; does not mean a construct functions in SNIPR."}
        for name, value in (("candidates", candidates), ("protein_diagnostics", diagnostics), ("construct_audit", list(constructs.values())), ("observations", observations), ("summary", summary), ("conflicts", conflicts), ("evidence_links", links)):
            write_json(run_dir / (name + ".json"), value)
        tsv(run_dir / "candidate_plan.tsv", candidates, CANDIDATE_COLUMNS)
        criteria = [dict(candidate_id=c["candidate_id"], protein_id=c["protein_id"], **r) for c in candidates for r in c.get("criteria_review", [])]
        tsv(run_dir / "criteria_review.tsv", criteria, ["candidate_id", "protein_id", "kind", "status", "interpretation", "records"])
        tsv(run_dir / "assay_observations.tsv", observations, ["observation_id", "construct_id", "endpoint", "raw_value", "unit", "normalized_value", "normalized_unit", "eligibility", "issues", "context", "source_ref"])
        tsv(run_dir / "construct_audit.tsv", list(constructs.values()), ["construct_id", "protein_id", "scaffold_id", "scaffold_version", "start", "end", "audit_status", "audit_issues"])
        with (run_dir / "candidate_fragments.fasta").open("w") as handle:
            for c in candidates:
                if c["sequence"] and c["design_status"] != "blocked":
                    handle.write(f">{c['candidate_id']} {c['protein_id']}:{c['start']}-{c['end']} UNVALIDATED\n{c['sequence']}\n")
        with (run_dir / "reviewed_fusions.fasta").open("w") as handle:
            for c in candidates:
                if c["fusion_sequence"]:
                    handle.write(f">{c['candidate_id']} UNVALIDATED synthetic={c['synthetic_only']}\n{c['fusion_sequence']}\n")
        for name, values in (("claims.jsonl", claims), ("evidence_records_v2.jsonl", project.get("evidence_records_v2", []))):
            (run_dir / name).write_text("".join(canonical(v) + "\n" for v in values), encoding="utf-8")
        (run_dir / "report.md").write_text(report(summary, candidates, diagnostics), encoding="utf-8")
        state["state"] = "complete"
        event("export", "complete")
        files = {p.name: file_hash(p) for p in sorted(run_dir.iterdir()) if p.is_file() and p.name != "manifest.json"}
        write_json(run_dir / "manifest.json", dict(state, artifacts=files, finished_at=now()))
        return {"run_dir": str(run_dir.resolve()), "execution": "computed", "summary": summary}
    except Exception as exc:
        state["state"] = "failed"
        event("failure", "error", error=str(exc))
        raise


def report(summary, candidates, diagnostics):
    lines = ["# Antigen-Receiver candidate review", "", "Configuration: Antibody-Sender -> Antigen-Receiver.", "", 
             f"Proteins: {summary['proteins']}; candidates: {summary['candidates']}; reviewed fusion sequences: {summary['fusion_sequences']}.",
             "", "These are proposed constructs, not a functional success ranking. Empty fusion FASTA is intentional when scaffold/review is missing.", "",
             "## Decisions", "", "| Protein | Candidate | Range | Design | Assembly | Review reasons |", "|---|---|---|---|---|---|"]
    for c in candidates:
        codes = ", ".join(sorted({r["code"] for r in c["risks"] + c["scaffold_issues"]}))
        lines.append(f"| {c['protein_id']} | {c['candidate_id']} | {c.get('start')}-{c.get('end')} | {c['design_status']} | {c['assembly_status']} | {codes} |")
    lines += ["", "## No supported candidate", ""]
    for d in diagnostics:
        if not d["candidate_count"]:
            lines.append(f"- {d['protein_id']}: " + ", ".join(r["code"] for r in d["diagnostics"]))
    lines += ["", "## Four independent endpoints", "", "Surface expression; recognition retention; basal activation; induced response.",
              "No automatic pooling or conversion of SNIPR Assay Y/N to a success label. See assay_observations.tsv for unresolved definitions and contexts.",
              "", "## Next experimental decisions", "", "Review full ECD versus truncated epitope scope; confirm actual scaffold and junctions; associate each construct with its sender and assay conditions. Validate all four endpoints independently.",
              "", "## Evidence boundaries", "", "Sequence-retained epitope is not binding-retained epitope. Potential glycosylation is a sequence motif, not observed occupancy. Length and cysteine flags are configurable review heuristics, not validated performance thresholds.",
              "", "See criteria_review.tsv for processing, shedding, junction cleavage, oligomerization, aggregation and culture-interference evidence. Missing means not assessed, not low risk. Reported native biology does not predict receiver background. Full antigen forms are preferred; optional domain fragments may lose screening epitopes.",
              "", "Private input snapshots and outputs must remain outside Git and skill installation bundles.", ""]
    return "\n".join(lines)
