"""Endpoint-specific receiver evidence, independent of fragment design classes.

No gene lists, fitted cutoffs, risk-count scores or automatic outcome labels.
The absence of an annotation signal never constitutes a negative assay result.
"""

from collections import Counter
from pathlib import Path

from .common import digest, file_hash, now, read_json, tsv, write_json
from .observations import ENDPOINTS, normalize_observation

POLICY_VERSION = "1.0"
UNKNOWN = "undetermined_not_low_risk"
LIMIT = "Design feasibility does not establish receiver expression, recognition, low background or inducibility."

# These are review questions, not a causal mapping from features to failures.
AXES = (
    ("folding_trafficking", ("domain_cut", "disulfide_partner_removed", "cysteine_rich_folding_review", "free_thiol_odd_cysteine", "dense_glycosylation"), (),
     "Inspect retained structural units and surface expression; motifs do not establish folding or glycan occupancy."),
    ("recognition_repertoire", ("known_epitope_loss", "domain_fragment_independence_unproven"), (),
     "Measure recognition with independent binders; intact sequence does not prove retained epitopes."),
    ("attachment_geometry", ("type_ii_attachment_orientation_change", "gpi_anchor_replaced", "long_fragment_geometry_review", "short_fragment_geometry_review"), (),
     "Review fragment orientation and spacing in the actual receiver; no universal length cutoff predicts signaling."),
    ("processing_junctions", ("processed_chain_segments_spanned", "extracellular_processed_product_not_covered"), ("processing", "native_shedding", "junction_cleavage"),
     "Inspect the actual antigen-junction-regulatory-core context; native shedding alone does not predict basal activation."),
    ("assembly_interactions", ("native_heteromer_context_requires_review", "multichain_partner_required"), ("oligomerization", "aggregation"),
     "Review interaction context without equating native complex membership with obligatory ECD partner dependence or receptor clustering."),
    ("expression_and_assay_context", (), ("culture_interference",),
     "Record backbone, host, expression regime, reporter, sender, timepoint and controls; do not transfer between receptor variants."),
)

ENDPOINT_ACTIONS = {
    "surface_expression": "Measure surface-localized receptor with an explicit live/singlet/receiver denominator; total expression is not surface expression.",
    "recognition_retention": "Measure specific binding to the actual receiver, with nonbinding controls; separate binding from reporter activation.",
    "basal_activation": "Measure unstimulated and nonbinding-control background with matched receptor expression and reporter-only/background controls.",
    "induced_response": "Measure binding-sender response and matched baseline; retain both raw values before calculating any induction contrast.",
}

FUNCTION_COLUMNS = ["candidate_id", "protein_id", "gene", "accession", "isoform", "is_primary",
                    "design_recommendation", "design_status", "assembly_status", "scaffold_id", "scaffold_version",
                    "functional_risk", "endpoint", "evidence_state", "measurement_count", "measurement_contexts",
                    "exact_measurements", "related_or_unresolved_records", "missing_information", "next_action"]


def _linked(candidate, construct):
    if not candidate.get("sequence"):
        return None
    antigen = (candidate.get("protein_id"), candidate.get("start"), candidate.get("end"), candidate.get("sequence"))
    recorded = (construct.get("protein_id"), construct.get("start"), construct.get("end"), construct.get("antigen_sequence"))
    if antigen != recorded:
        return None
    exact = (candidate.get("fusion_sequence") and candidate.get("scaffold_id") and candidate.get("scaffold_version")
             and candidate["fusion_sequence"] == construct.get("fusion_sequence")
             and candidate["scaffold_id"] == construct.get("scaffold_id")
             and candidate["scaffold_version"] == construct.get("scaffold_version")
             and construct.get("audit_status") == "sequence_consistent_not_functionally_validated")
    return "full_fusion_exact" if exact else "antigen_only_not_equivalent_receptor"


def _measurement_context(record):
    return {k: record.get(k) for k in ("construct_id", "endpoint", "context", "normalized_unit", "gate_path", "channel", "denominator", "statistic")}


def assess_receiver(candidate, observations=(), constructs=None):
    """Associate existing measurements, never turn them into inferred success labels."""
    constructs = constructs or {}
    linked = {cid: _linked(candidate, c) for cid, c in constructs.items()}
    risks = candidate.get("risks", [])
    axes = []
    for name, codes, kinds, question in AXES:
        signals = [{"origin": "candidate_design_annotation_or_heuristic", "record": r} for r in risks if r.get("code") in codes]
        context_records = []
        for criterion in candidate.get("criteria_review", []):
            if criterion.get("kind") not in kinds:
                continue
            for record in criterion.get("records", []):
                context_records.append(record)
                # Do not promote missing, unresolved or out-of-fragment reports.
                if record.get("effective_status") == "reported" and not record.get("interpretation", "").startswith("reported_region_not_retained"):
                    signals.append({"origin": "sourced_context_report_not_function_prediction", "record": record})
        axes.append({"axis": name, "state": "review_signal_present" if signals else "not_established",
                     "signals": signals, "context_records": context_records, "next_action": question,
                     "functional_prediction": "not_made"})
    endpoints = {}
    for endpoint in ENDPOINTS:
        exact, related = [], []
        for observation in observations:
            if observation.get("endpoint") != endpoint or not linked.get(observation.get("construct_id")):
                continue
            link = linked[observation["construct_id"]]
            if link == "full_fusion_exact" and observation.get("eligibility") == "eligible_measurement_not_success_label":
                exact.append(dict(observation))
            else:
                related.append({"link_scope": link, "observation": observation})
        groups = {}
        for record in exact:
            key = digest(_measurement_context(record))
            group = groups.setdefault(key, {"context": _measurement_context(record), "values": [], "observation_ids": []})
            group["values"].append(record.get("normalized_value"))
            group["observation_ids"].append(record.get("observation_id"))
        for group in groups.values():
            group["conflict"] = len(set(group["values"])) > 1
        conflicting = any(g["conflict"] for g in groups.values())
        state = ("conflicting_measurements" if conflicting else "measured_in_recorded_context") if exact else \
            "related_or_unresolved_only" if related else "not_measured"
        missing = ["predeclared_endpoint_acceptance_criteria_and_independent_validation"]
        if not candidate.get("fusion_sequence"):
            missing.append("confirmed_complete_receiver_sequence_and_backbone")
        if not exact:
            missing.append("eligible_same_receiver_endpoint_measurement")
        if endpoint in {"basal_activation", "induced_response"}:
            missing.append("matched_expression_and_control_comparison_not_computed")
        if candidate.get("synthetic_only"):
            missing.append("synthetic_fixture_not_biological_validation")
        endpoints[endpoint] = {"evidence_state": state, "functional_risk": UNKNOWN,
            "measurement_count": len(exact), "measurement_contexts": list(groups.values()),
            "exact_measurements": exact, "related_or_unresolved_records": related,
            "missing_information": missing, "next_action": ENDPOINT_ACTIONS[endpoint]}
    return {"policy_version": POLICY_VERSION, "scope": "candidate_x_receiver_x_assay_context",
            "functional_risk": UNKNOWN, "functional_compatibility": "not_established",
            "mechanism_review": axes, "endpoints": endpoints,
            "generic_warnings_not_failure_predictions": [r for r in risks if r.get("code") in {
                "contextual_risk_evidence_missing", "epitope_coverage_unknown", "prediction_requires_annotation_review", "potential_n_glycosylation"}],
            "claim_limit": LIMIT, "success_probability": None, "calibration": "not_performed",
            "anti_leakage_policy": "No gene-specific failure rules; observed cases used for rule development cannot be a held-out test set."}


def function_summary(candidates):
    usable = [c for c in candidates if c.get("sequence") and c.get("design_status") != "blocked"]
    return {"policy_version": POLICY_VERSION, "candidate_records": len(candidates),
            "usable_candidate_records": len(usable), "denominator_unit": "candidate_not_gene_or_experiment",
            "risk_states": dict(Counter(c["receiver_function"]["functional_risk"] for c in usable)),
            "endpoint_evidence_states": {e: dict(Counter(c["receiver_function"]["endpoints"][e]["evidence_state"] for c in usable)) for e in ENDPOINTS},
            "verified_functional_compatibility": "not_assessed", "claim_limit": LIMIT}


def export_receiver_function(outdir, candidates):
    rows, mechanisms = [], []
    for c in candidates:
        assessment = c["receiver_function"]
        metadata = {k: c.get(k) for k in ("candidate_id", "protein_id", "gene", "accession", "isoform", "is_primary", "design_status", "assembly_status", "scaffold_id", "scaffold_version")}
        metadata["design_recommendation"] = c.get("screening_recommendation")
        for endpoint, value in assessment["endpoints"].items():
            rows.append(dict(metadata, endpoint=endpoint, **value))
        for axis in assessment["mechanism_review"]:
            mechanisms.append(dict(candidate_id=c["candidate_id"], **axis))
    outdir = Path(outdir)
    tsv(outdir / "receiver_function.tsv", rows, FUNCTION_COLUMNS)
    tsv(outdir / "receiver_mechanism_review.tsv", mechanisms,
        ["candidate_id", "axis", "state", "signals", "context_records", "functional_prediction", "next_action"])
    summary = function_summary(candidates)
    write_json(outdir / "receiver_function_summary.json", summary)
    (outdir / "RECEIVER_FUNCTION.md").write_text(
        "# 片段设计与受体功能必须分开\n\n"
        "常规候选 = 能提出有依据的片段，不等于低自激活、表达正常或可诱导。条件性候选也不等于功能失败。\n\n"
        "receiver_function.tsv 每个候选分四类终点展示。undetermined_not_low_risk 表示尚不能判定功能风险，绝非低风险。"
        "measured_in_recorded_context 只表示同完整构建有可解释测量，不是通过验收或跨条件有效。\n\n"
        "receiver_mechanism_review.tsv 中的结构/方向/加工/相互作用线索只是待核验问题，不能计作失败预测命中；"
        "没有线索仍是 not_established，不是阴性。宿主、表达制度、骨架、连接和对照不能由原蛋白注释替代。\n\n"
        "无成功率、无未经验证的高低风险分级、无基因黑名单、无失败表反向调阈值。"
        "失败案例仅作开发问题来源；未来验证必须冻结规则，按蛋白/家族及完整构建分组，保留独立成功与失败对照。\n\n"
        f"候选记录：{summary['candidate_records']}；可用片段记录：{summary['usable_candidate_records']}。"
        "计数单位是候选，不是基因或独立实验。具体来源、上下文与冲突保留在 JSON/TSV。\n", encoding="utf-8")
    return summary


def audit_run(run_dir, outdir, resume=False):
    """Add an immutable overlay to a verified historical bundle; never rescreen it."""
    from . import __version__
    from .harness import engine_hash, verify_bundle
    source = Path(run_dir).resolve()
    if not verify_bundle(source):
        raise ValueError("Source run integrity check failed")
    root = Path(outdir).resolve()
    if root == source or source in root.parents:
        raise ValueError("Receiver audit output must be outside the frozen source bundle")
    code_hash = engine_hash()
    input_hash = file_hash(source / "manifest.json")
    key = digest({"source_manifest_sha256": input_hash, "engine_sha256": code_hash, "policy": POLICY_VERSION})
    dest = root / key[:20]
    if dest.exists():
        if resume and verify_bundle(dest):
            return {"run_dir": str(dest), "execution": "verified_cache_hit", "summary": read_json(dest / "receiver_function_summary.json")}
        raise FileExistsError("Audit exists or is corrupt; preserved. Choose a new output root.")
    candidates = read_json(source / "candidates.json")
    constructs = {c["construct_id"]: c for c in read_json(source / "construct_audit.json")} if (source / "construct_audit.json").exists() else {}
    raw = read_json(source / "observations.json") if (source / "observations.json").exists() else []
    observations = [normalize_observation(o, constructs) for o in raw]
    for candidate in candidates:
        candidate["receiver_function"] = assess_receiver(candidate, observations, constructs)
    if not verify_bundle(source) or file_hash(source / "manifest.json") != input_hash:
        raise ValueError("Source changed during audit")
    dest.mkdir(parents=True)
    write_json(dest / "candidates.json", candidates)
    summary = export_receiver_function(dest, candidates)
    write_json(dest / "source_link.json", {"source_run": str(source), "source_manifest_sha256": input_hash,
        "source_engine_sha256": read_json(source / "manifest.json").get("engine_sha256"),
        "design_rules_changed": False, "rescreened": False, "additional_workbooks_read": False,
        "data_privacy": "inherits_source_bundle_never_assumed_public"})
    artifacts = {p.name: file_hash(p) for p in sorted(dest.iterdir()) if p.is_file()}
    write_json(dest / "manifest.json", {"state": "complete", "version": __version__, "engine_sha256": code_hash,
        "run_key": key, "source_manifest_sha256": input_hash, "finished_at": now(), "artifacts": artifacts})
    return {"run_dir": str(dest), "execution": "computed", "summary": summary}
