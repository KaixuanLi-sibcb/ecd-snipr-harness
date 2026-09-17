"""Native claim ledger plus lossless legacy-v2 interoperability."""

from .common import digest


def claims_for(project, candidates, protein_diagnostics, observations):
    claims = []
    for p in project["proteins"]:
        claims.append({"claim_type": "reference_sequence", "subject_id": p["protein_id"],
                       "status": "recorded" if p.get("sequence") else "missing", "evidence_class": p.get("evidence", {}).get("kind", "unknown"),
                       "source": p.get("evidence", {}), "value_sha256": digest(p.get("sequence", "")),
                       "limitations": "Recorded source is not proof of expression or SNIPR function"})
    for c in candidates:
        claims.append({"claim_type": "candidate_design", "subject_id": c["candidate_id"], "protein_id": c["protein_id"],
                       "status": "computed" if c["sequence"] else "error", "evidence_class": "deterministic_rule",
                       "source": c.get("boundary_evidence", {}), "review_key": c["review_key"],
                       "value": {"start": c.get("start"), "end": c.get("end"), "sequence_sha256": digest(c["sequence"]), "risks": c["risks"]}})
        for criterion in c.get("criteria_review", []):
            claims.append({"claim_type": "contextual_criterion", "subject_id": c["candidate_id"], "status": "missing" if criterion["status"] == "missing" else "recorded", "evidence_class": "scoped_source_records", "value": criterion, "limitations": "No cross-context transfer, absence-to-negative conversion or functional probability"})
    for d in protein_diagnostics:
        if not d["candidate_count"]:
            claims.append({"claim_type": "candidate_availability", "subject_id": d["protein_id"], "status": "missing", "evidence_class": "deterministic_rule", "reasons": d["diagnostics"]})
    for o in observations:
        claims.append({"claim_type": "assay_observation", "subject_id": o["observation_id"], "construct_id": o.get("construct_id", ""),
                       "status": "recorded" if o["eligibility"] == "eligible_measurement_not_success_label" else "missing" if o["eligibility"] == "missing" else "unresolved",
                       "evidence_class": "local_experiment", "context": o.get("context", {}), "endpoint": o.get("endpoint", "unconfirmed"),
                       "source": o.get("source_ref", {}), "value": o.get("normalized_value"), "issues": o["issues"]})
    for claim in claims:
        claim["claim_id"] = digest(claim)
    links = [{"record_sha256": digest(r), "legacy_candidate_id": r["candidate_id"], "evidence_type": r["evidence_type"],
              "link_status": "unresolved_unless_explicitly_mapped", "note": "Legacy gene/candidate identity alone does not resolve construct or assay context"} for r in project.get("evidence_records_v2", [])]
    conflicts = []
    groups = {}
    for o in observations:
        # Replicate, batch, stimulus and sender remain part of the comparison key.
        key = digest([o.get("construct_id"), o.get("endpoint"), o.get("context"), o.get("unit"), o.get("statistic"), o.get("channel"), o.get("gate_path"), o.get("denominator")])
        groups.setdefault(key, []).append(o)
    for key, values in groups.items():
        nonmissing = [o for o in values if o.get("raw_value") is not None and str(o["raw_value"]).strip()]
        if len({str(o["raw_value"]) for o in nonmissing}) > 1:
            conflicts.append({"context_key": key, "type": "duplicate_context_value_conflict", "observations": [o["observation_id"] for o in nonmissing], "resolution": "pending_no_overwrite"})
    return claims, links, conflicts
