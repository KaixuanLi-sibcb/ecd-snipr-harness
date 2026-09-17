"""Keep construct/sender/condition identities and four assay endpoints separate."""

import math
from .common import digest

ENDPOINTS = ("surface_expression", "recognition_retention", "basal_activation", "induced_response")


def normalize_observation(record, constructs):
    result = dict(record)
    result.update(normalized_value=None, eligibility="unresolved", issues=[])
    issues = result["issues"]
    result["observation_key"] = digest(record)
    value = record.get("raw_value")
    if value is None or str(value).strip().lower() in {"", "na", "n/a", "x", "pending"}:
        result["eligibility"] = "missing"
        issues.append("blank_or_placeholder_is_not_negative")
        return result
    if record.get("endpoint") not in ENDPOINTS:
        issues.append("endpoint_definition_unconfirmed")
    if record.get("construct_id") not in constructs:
        issues.append("construct_identity_unresolved")
    elif constructs[record["construct_id"]]["audit_status"] != "sequence_consistent_not_functionally_validated":
        issues.append("construct_sequence_or_scaffold_requires_review")
    context = record.get("context", {})
    required = ("batch_id", "host_cell", "sender_id", "antibody_id", "stimulus", "condition_id", "timepoint", "replicate_id")
    for k in required:
        if not context.get(k) or context.get(k) == "unknown":
            issues.append("context_missing:" + k)
    if record.get("endpoint") == "basal_activation":
        if context.get("stimulus") not in {"none", "nonbinding_control"}:
            issues.append("basal_condition_not_confirmed")
    if record.get("endpoint") in {"induced_response", "recognition_retention"}:
        if context.get("stimulus") != "binding_sender":
            issues.append("stimulated_condition_not_confirmed")
    if not record.get("source_ref"):
        issues.append("source_cell_or_record_missing")
    if not record.get("gate_path") or not record.get("channel"):
        issues.append("gating_or_channel_unconfirmed")
    try:
        if isinstance(value, bool):
            raise ValueError()
        number = float(str(value).strip())
        if not math.isfinite(number):
            raise ValueError()
    except (ValueError, TypeError):
        issues.append("not_a_numeric_measurement_no_Y_N_conversion")
        return result
    unit = record.get("unit")
    if unit in {"percent", "fraction"}:
        if not record.get("denominator"):
            issues.append("fraction_denominator_unconfirmed")
        limit = 100 if unit == "percent" else 1
        if not 0 <= number <= limit:
            issues.append("out_of_range_for_declared_unit")
        else:
            result["normalized_value"] = number / limit
            result["normalized_unit"] = "fraction"
    elif unit == "fluorescence_au":
        if record.get("statistic") not in {"mean", "median", "geometric_mean"}:
            issues.append("MFI_statistic_unconfirmed")
        if number < 0:
            issues.append("negative_value_requires_instrument_review")
        result["normalized_value"] = number
        result["normalized_unit"] = unit
    else:
        issues.append("unit_unconfirmed_no_automatic_guessing")
    result["eligibility"] = "eligible_measurement_not_success_label" if not issues else "unresolved"
    return result


def assay_summary(observations, constructs):
    eligible = [o for o in observations if o["eligibility"] == "eligible_measurement_not_success_label"]
    return {
        "records": len(observations), "eligible_measurements": len(eligible),
        "distinct_construct_ids": len(constructs),
        "sequence_consistent_constructs": sum(c["audit_status"] == "sequence_consistent_not_functionally_validated" for c in constructs.values()),
        "unique_full_fusion_sequences": len({c["fusion_sequence"] for c in constructs.values() if c.get("fusion_sequence")}),
        "eligible_constructs_by_endpoint": {e: len({o["construct_id"] for o in eligible if o["endpoint"] == e}) for e in ENDPOINTS},
        "batches": sorted({o.get("context", {}).get("batch_id", "unknown") for o in observations}),
        "success_failure_class_distribution": "not_defined",
        "model_readiness": "not_assessed_no_functional_classifier_or_probability",
        "warning": "No pooling across senders, constructs, conditions or batches; no causal inference from coverage counts.",
    }
