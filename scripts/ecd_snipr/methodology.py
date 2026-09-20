"""Orthogonal evidence, repertoire tradeoffs and receiver review planning.

These records describe annotation support and design choices, never calibrated
confidence or functional probabilities. Unknown ECO codes remain unclassified.
"""

from collections import Counter
from .common import digest

METHOD_VERSION = "1.0"
ECO_CLASSES = {
    "ECO:0000269": "experimental_annotation",
    "ECO:0000305": "curator_inference",
    "ECO:0000250": "similarity_transfer",
    "ECO:0000255": "sequence_model_inference",
    "ECO:0000256": "automatic_assertion",
    "ECO:0000259": "sequence_model_inference",
    "ECO:0000303": "author_statement",
    "ECO:0000313": "imported_information",
}
FORM_PRIORITY = {"full_ecd": 0, "mature_gpi": 1, "mature_secreted": 2,
                 "shed_product": 3, "domain_fragment": 4}
SELECTION_AXES = ["annotated_integrity_disruption", "mapped_epitope_loss",
                  "antigen_form_priority", "candidate_id_tiebreak"]


def evidence_support(evidence):
    evidence = evidence or {}
    raw = evidence.get("eco") or []
    if isinstance(raw, (str, dict)):
        raw = [raw]
    codes = sorted({e.get("evidenceCode", "") if isinstance(e, dict) else str(e)
                    for e in raw if e})
    codes = [c for c in codes if c]
    kind = evidence.get("kind", "missing")
    if kind == "synthetic_fixture":
        labels = ["synthetic_only"]
    elif codes:
        labels = sorted({ECO_CLASSES.get(c, "unclassified_eco") for c in codes})
    else:
        labels = [{"curated_annotation": "curated_support_unspecified",
                   "prediction": "prediction_method_unspecified",
                   "local_experiment": "local_experiment_context_requires_review"}.get(kind, "missing")]
    return {"support_types": labels, "eco_codes": codes, "provenance": evidence,
            "interpretation": "Annotation-level provenance, not SNIPR evidence or a confidence score"}


def _feature_record(f):
    return {"kind": f.get("kind"), "name": f.get("name", ""),
            "start": f.get("start"), "end": f.get("end"),
            "support": evidence_support(f.get("evidence"))}


def annotation_support(protein, candidate):
    fs = protein.get("features", [])
    start, end = candidate.get("start"), candidate.get("end")
    valid = type(start) is int and type(end) is int
    relevant_domains = [f for f in fs if f.get("kind") == "domain" and valid
                        and type(f.get("start")) is int and type(f.get("end")) is int
                        and f["start"] <= end and start <= f["end"]]
    return {
        "reference": evidence_support(protein.get("evidence")),
        "boundary": evidence_support(candidate.get("boundary_evidence")),
        "topology": [_feature_record(f) for f in fs if f.get("kind") in
                     {"extracellular", "transmembrane", "gpi_attachment_site", "signal_peptide", "chain"}],
        "overlapping_domains": [_feature_record(f) for f in relevant_domains],
        "domain_coverage_status": "annotations_present" if relevant_domains else "not_assessed_no_overlapping_annotation",
        "reviewed_entry_is_experimental_proof": False,
        "note": "No evidence propagation between components; absent ECO is not experimental support",
    }


def candidate_tradeoff(candidate):
    risks = {r["code"] for r in candidate.get("risks", [])}
    ep = candidate.get("epitope_review", [])
    loss = [e for e in ep if e.get("state") != "sequence_retained_not_binding_proven"]
    integrity = sorted(risks & {"domain_cut", "disulfide_partner_removed"})
    return {
        "integrity_disruptions": integrity,
        "mapped_epitope_status": "loss_or_partial_loss" if loss else "sequence_retained_not_binding_proven" if ep else "unknown",
        "mapped_epitope_losses": loss,
        "domains_retained": candidate.get("domains_retained", []),
        "domains_cut": candidate.get("domains_cut", []),
        "domains_omitted": candidate.get("domains_omitted", []),
        "disulfide_crossings": (candidate.get("molecular_profile") or {}).get("disulfides_partial", []),
        "repertoire_statement": "Omitted domains/epitopes narrow the tested antigen repertoire; unknown epitopes remain unknown",
        "selection_axes": SELECTION_AXES,
        "selection_key": [bool(integrity), bool(loss), FORM_PRIORITY.get(candidate.get("antigen_form_type"), 9),
                          candidate["candidate_id"]],
        "selection_limit": "Explicit heuristic, not learned/calibrated. Unknown epitope coverage is not proven retention",
    }


def selection_key(candidate):
    return tuple(candidate_tradeoff(candidate)["selection_key"])


def receiver_review(candidate):
    """Candidate-context questions, not an assumed laboratory scaffold."""
    codes = {r["code"] for r in candidate.get("risks", [])}
    topology = candidate.get("topology", "unknown")
    checks = [
        ("receiver_architecture", "not_evaluated", "Provide the actual antigen slot, regulatory modules, junctions, SP, host and reporter; do not substitute sender PDGFR-LC", []),
        ("attachment_geometry", "requires_review", "If the antigen slot is N-terminal to receiver TM, attachment is at the fragment C terminus; verify native orientation and steric accessibility", ["type_ii_attachment_orientation_change"] if topology == "type_ii" else []),
        ("recognition_repertoire", "requires_review", "Compare retained/omitted domains and mapped epitopes; sequence retention does not establish binding or new-antibody repertoire", sorted(codes & {"known_epitope_loss", "domain_cut", "domain_omitted_epitope_scope_changed", "epitope_coverage_unknown"})),
        ("folding_and_partners", "requires_review", "Review domain boundaries, disulfide crossings and native complex context without assuming obligatory partner dependence", sorted(codes & {"disulfide_partner_removed", "domain_cut", "native_heteromer_context_requires_review"})),
        ("processing_and_junctions", "not_evaluated", "Native processing/shedding annotations do not establish cleavage of a new fusion; assess actual junction sequence and cell context", sorted(codes & {"gpi_anchor_replaced", "multiple_processed_chains_dependency_unknown", "shed_product_tethering_unvalidated"})),
    ]
    return {"configuration": "antibody_sender_antigen_receiver",
            "status": "planning_only_scaffold_unverified",
            "candidate_design_status": candidate.get("design_status"),
            "checks": [{"check": k, "status": s, "question": q, "candidate_reason_codes": c} for k, s, q, c in checks],
            "endpoints": {k: "not_measured" for k in ("surface_expression", "recognition_retention", "basal_activation", "induced_response")},
            "required_observation_key": ["construct_id", "candidate_id", "scaffold_version", "sender_or_antibody_id", "batch", "replicate", "condition", "gating_denominator", "unit"],
            "fusion_sequence_authorized": False}


def route_diagnostic(record, candidates):
    cls = record.get("screening_recommendation")
    codes = set(record.get("reason_codes", []))
    if record.get("processing_status") == "technical_failure":
        state, action = "technical_failure", "Repair fetch/cache/input integrity, then retry this record"
    elif record.get("duplicate_of"):
        state, action = "duplicate_input", "Consult the first occurrence; do not count another reference"
    elif record.get("scope_status") == "out_of_scope":
        state, action = "outside_declared_scope", "Separate natural-surface scope from exploratory engineered presentation; no impossibility claim"
    elif cls == "insufficient_evidence":
        state, action = "core_annotation_gap_or_conflict", "Resolve identity, reference sequence, topology or mature boundaries indicated by reason_codes"
    elif cls == "no_standard_route":
        if any(c.startswith("lab_rule_") for c in codes) or any(r.get("effect") == "downgrade_to_no_standard_route" for r in record.get("lab_rules", {}).get("applied", [])):
            state, action = "lab_context_route_restriction", "Review the sourced laboratory rule in its original construct/assay context"
        elif candidates and all(c.get("design_status") == "blocked" for c in candidates):
            state, action = "proposed_routes_blocked", "Review exact blocking reasons; do not bypass coordinate or native-sequence exclusions"
        else:
            state, action = "no_supported_route_in_current_annotations", "Seek sourced autonomous extracellular domains/mature forms or another presentation route; absent annotation is not negative evidence"
    elif cls in {"standard_candidate", "conditional_candidate"}:
        state, action = "supported_candidate_unvalidated", "Review candidate comparison and receiver checklist before actual scaffold assembly and four-endpoint experiments"
    else:
        state, action = "not_assessed", "Inspect unresolved processing/scope status"
    return {"state": state, "reason_codes": sorted(codes), "next_action": action,
            "biological_impossibility_claim": False, "method_version": METHOD_VERSION}


def review_queue(records, candidates, per_stratum=2):
    """Deterministic purposive audit sample. It cannot estimate accuracy."""
    if type(per_stratum) is not int or not 1 <= per_stratum <= 100:
        raise ValueError("review-per-stratum must be an integer between 1 and 100")
    by_id = {c["candidate_id"]: c for c in candidates}
    grouped = {}
    for r in records:
        if r.get("duplicate_of"):
            continue
        c = by_id.get(r.get("primary_candidate_id"), {})
        support = (c.get("annotation_support", {}).get("boundary", {}).get("support_types", ["not_available"]))
        tradeoff = c.get("candidate_comparison", {})
        stratum = [r.get("scope_status"), r.get("topology"), r.get("screening_recommendation"),
                   r.get("route_diagnostic", {}).get("state"), support, tradeoff.get("integrity_disruptions", []),
                   tradeoff.get("mapped_epitope_status", "not_available")]
        grouped.setdefault(digest(stratum), {"stratum": stratum, "records": []})["records"].append(r)
    rows = []
    for key, group in sorted(grouped.items()):
        selected = sorted(group["records"], key=lambda r: digest([r.get("accession"), r.get("isoform"), r["entry_id"]]))[:per_stratum]
        for r in selected:
            rows.append({"entry_id": r["entry_id"], "accession": r.get("accession", ""), "isoform": r.get("isoform", ""),
                         "candidate_id": r.get("primary_candidate_id", ""), "stratum_id": key[:16],
                         "stratum": group["stratum"], "stratum_population": len(group["records"]),
                         "selection_policy": "deterministic_purposive_not_accuracy_sample",
                         "review_status": "pending", "independent_boundary": "", "reviewer": "", "source": "",
                         "disagreement_reason": "", "surface_expression": "", "recognition_retention": "",
                         "basal_activation": "", "induced_response": ""})
    return rows


def methodology_summary(records, candidates):
    primary = [c for c in candidates if c.get("is_primary")]
    return {"version": METHOD_VERSION,
            "route_states": dict(Counter(r["route_diagnostic"]["state"] for r in records)),
            "route_denominator": len(records), "route_denominator_definition": "All preserved input rows; duplicates explicitly separate",
            "primary_boundary_support_combinations": dict(Counter("+".join(c["annotation_support"]["boundary"]["support_types"]) for c in primary)),
            "primary_denominator": len(primary),
            "primary_with_integrity_disruption": sum(bool(c["candidate_comparison"]["integrity_disruptions"]) for c in primary),
            "primary_with_mapped_epitope_loss": sum(bool(c["candidate_comparison"]["mapped_epitope_losses"]) for c in primary),
            "primary_with_unknown_epitopes": sum(c["candidate_comparison"]["mapped_epitope_status"] == "unknown" for c in primary),
            "validation_status": "computational_only_not_biologically_validated"}
