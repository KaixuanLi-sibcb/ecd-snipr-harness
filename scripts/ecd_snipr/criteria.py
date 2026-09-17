"""Contextual evidence review, deliberately not a receptor success classifier."""

from .common import interval, overlaps, sourced

CONSIDERATIONS = (
    "processing", "native_shedding", "junction_cleavage", "oligomerization",
    "aggregation", "culture_interference",
)
FORMS = {"full_ecd", "mature_secreted", "mature_gpi", "shed_product", "domain_fragment"}
STATUSES = {"reported", "not_observed_in_context", "missing", "not_applicable"}
SCOPES = {"native_protein", "receiver_construct", "assay_context"}


def assess_context(protein, bounds):
    """Retain every report/conflict and its scope; absence is never low risk."""
    rows, risks = [], []
    for kind in CONSIDERATIONS:
        records = [r for r in protein.get("considerations", []) if r["kind"] == kind]
        if not records:
            rows.append({"kind": kind, "status": "missing", "interpretation": "not_assessed", "records": []})
            continue
        states, evaluated = [], []
        for record in records:
            state = record["status"]
            interpretation = "requires_context_review"
            if state == "missing":
                interpretation = "not_assessed"
            elif not sourced(record) or not record.get("reason") or not record.get("context"):
                state, interpretation = "unresolved", "source_reason_or_context_missing"
            if "start" in record or "end" in record:
                try:
                    region = interval(record, len(protein["sequence"]))
                    if not overlaps(bounds, region):
                        interpretation = "reported_region_not_retained_not_proof_of_no_effect"
                except ValueError:
                    state, interpretation = "unresolved", "invalid_context_interval"
            if record.get("scope") != "native_protein" and not all(record.get("context", {}).get(k) for k in ("construct_id", "condition_id")):
                state, interpretation = "unresolved", "construct_or_condition_identity_missing"
            if state == "not_observed_in_context":
                if record.get("evidence", {}).get("kind") == "prediction":
                    state, interpretation = "unresolved", "prediction_is_not_a_negative_observation"
                else:
                    interpretation = "negative_observation_only_in_recorded_context_not_general_safety"
            if state == "not_applicable":
                interpretation = "user_scoped_exclusion_requires_review"
            if state == "reported" and kind in {"native_shedding", "processing", "junction_cleavage"}:
                interpretation += ";not_a_prediction_of_receiver_basal_activation"
            evaluated.append(dict(record=record, effective_status=state, interpretation=interpretation))
            states.append(state)
        distinct = set(states) - {"missing"}
        status = "mixed_context_records" if len(distinct) > 1 else next(iter(distinct), "missing")
        rows.append({"kind": kind, "status": status, "interpretation": "no_cross_context_pooling", "records": evaluated})
        if distinct:
            risks.append({"code": kind + "_context_review", "severity": "review", "detail": status})
    return rows, risks
