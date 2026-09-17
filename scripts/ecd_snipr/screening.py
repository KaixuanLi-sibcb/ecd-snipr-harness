"""Independent screening-recommendation layer and the batch screening runner.

`screening_recommendation` is separate from design_status / assembly_status /
functional_status. A simple transparent decision table assigns one of four
classes: standard_candidate, conditional_candidate, no_standard_route (a route
limitation, not "protein unusable forever") or insufficient_evidence (missing
core information, not a biological failure). Out-of-scope objects are recorded
via scope_status; technical failures via processing_status. No composite
scores: every recommendation carries reason_codes, a short rationale, the
evidence used and what is missing.
"""

import time
from collections import Counter
from pathlib import Path
from . import __version__
from .acquisition import disposition, load_set_proteins
from .common import digest, file_hash, now, read_json, tsv, write_json, canonical
from .design import propose
from .harness import CANDIDATE_COLUMNS, engine_hash, verify_bundle

CLASSES = ("standard_candidate", "conditional_candidate", "no_standard_route", "insufficient_evidence")

# Issues that make a supported candidate conditional rather than standard:
# a specific, named design question that needs special design or verification.
CONDITIONAL_CODES = {
    "type_ii_attachment_orientation_change",
    "gpi_anchor_replaced",
    "gpi_processing_boundary_unconfirmed",
    "secreted_antigen_scope_extension",
    "secreted_location_unconfirmed",
    "multiple_processed_chains_dependency_unknown",
    "discontinuous_or_multipass",
    "domain_fragment_independence_unproven",
    "prediction_requires_annotation_review",
    "plasma_membrane_not_confirmed",
}

# Core information missing or internally conflicting: screening cannot give a
# reliable recommendation. This is an evidence state, not a biological failure.
INSUFFICIENT_CODES = {
    "reference_sequence_invalid",
    "reference_provenance_missing",
    "isoform_unresolved",
    "topology_unknown",
    "feature_annotation_invalid",
    "mature_chain_boundary_missing",
    "single_pass_topology_incomplete_or_conflicting",
    "topology_orientation_conflict",
    "multipass_topology_incomplete",
    "secreted_topology_conflict",
}

FORM_PRIORITY = {"full_ecd": 0, "mature_gpi": 1, "mature_secreted": 2, "shed_product": 3, "domain_fragment": 4}

# Feature kinds whose coordinates define or block candidates. Fuzzy/invalid
# annotations of these kinds stay blocking when they are the only annotation
# of that kind (a genuine coordinate problem is preserved). Fuzzy secondary
# annotations (a shed-form chain on a type-I receptor, an alternative splice
# chain, a fuzzy domain/disulfide/epitope) are deferred with an explicit
# record instead of poisoning an otherwise exact candidate.
CANDIDATE_ESSENTIAL_KINDS = {"extracellular", "transmembrane", "signal_peptide", "cytoplasmic", "gpi_signal", "propeptide"}
CHAIN_ESSENTIAL_TOPOLOGIES = {"secreted", "gpi"}

CONTEXTUAL_KINDS = ("processing", "native_shedding", "junction_cleavage", "oligomerization", "aggregation", "culture_interference")

LAB_RULE_EFFECTS = {"annotate", "downgrade_to_conditional", "downgrade_to_no_standard_route"}

SCREENING_COLUMNS = [
    "entry_id", "protein_id", "gene", "accession", "isoform", "reference_basis",
    "experimental_isoform_confirmation", "input_row", "input_value", "duplicate_of",
    "membrane_set_membership", "natural_cell_surface_target", "scope_status", "scope_reason_codes",
    "extension_set", "topology", "location", "screening_recommendation",
    "primary_candidate_id", "primary_antigen_form", "primary_start", "primary_end", "primary_length",
    "reason_codes", "rationale", "main_risks", "missing_info", "not_evaluated", "deferred_annotations",
    "lab_rules_status", "processing_status", "processing_error",
]


def load_lab_rules(path):
    """Import sourced laboratory experience rules; never silent, never an upgrade."""
    rules = read_json(path)
    items = rules.get("rules") if isinstance(rules, dict) else rules
    if not isinstance(items, list):
        raise ValueError("Lab rules file must contain a rules list")
    for rule in items:
        if not rule.get("rule_id") or not isinstance(rule.get("match", {}), dict) or not rule["match"]:
            raise ValueError("Every lab rule needs rule_id and a non-empty match object")
        if set(rule["match"]) - {"gene", "accession", "topology"}:
            raise ValueError("Lab rules match on gene/accession/topology only")
        if rule.get("effect") not in LAB_RULE_EFFECTS:
            raise ValueError("Unknown lab rule effect: " + str(rule.get("effect")))
        evidence = rule.get("evidence", {})
        if evidence.get("kind") not in {"local_experiment", "curated_annotation"} or not evidence.get("source") or not evidence.get("version"):
            raise ValueError("Lab rule " + str(rule.get("rule_id")) + " lacks sourced provenance (kind/source/version)")
        if not rule.get("rationale"):
            raise ValueError("Lab rule " + str(rule.get("rule_id")) + " needs a rationale")
    return items


def _lab_rules_for(protein, rules):
    """Apply matching lab rules. Rules may annotate or downgrade, never upgrade."""
    applied = []
    for rule in rules or []:
        match = rule["match"]
        if any(protein.get(key) != value for key, value in match.items()):
            continue
        applied.append({"rule_id": rule["rule_id"], "effect": rule["effect"], "reason_code": rule.get("reason_code", ""),
                        "rationale": rule["rationale"], "evidence": rule["evidence"]})
    return applied


def junction_notes(candidate, protein):
    """What an experimenter must know to connect this fragment to the receiver.

    With no verified scaffold these are notes only; no fusion sequence is made.
    """
    notes = ["片段接入接收端抗原槽（antigen slot）；具体连接位置与 linker 须待真实骨架确认，当前不生成融合序列"]
    topology = protein.get("topology", "unknown")
    if topology == "type_ii":
        notes.append("II 型膜蛋白：天然为 N 端膜系留，接入本接收端后变为 C 端系留，接入方向改变需重点核验")
    if topology == "gpi":
        notes.append("GPI：候选保留 ω 残基，成熟加工边界与膜呈递需实验核对；ω+1 起的 GPI 信号不得保留")
    if candidate.get("antigen_form_type") == "domain_fragment":
        notes.append("结构域备选：完整结构域有来源，但独立折叠与表位覆盖需审阅")
    if candidate.get("antigen_form_type") == "mature_secreted":
        notes.append("分泌链扩展对象：系留呈递方式非常规天然状态，需专门设计")
    if candidate.get("origin") == "domain_alternative":
        notes.append("备选方案，非主要建议；主要建议见 primary_candidate_id")
    return notes


def _usable(candidate):
    return bool(candidate.get("sequence")) and candidate.get("design_status") != "blocked"


def pick_primary(candidates):
    """Explicit deterministic primary-candidate rule; alternates keep reasons."""
    usable = [c for c in candidates if _usable(c)]

    def key(c):
        codes = {r["code"] for r in c.get("risks", [])}
        return (FORM_PRIORITY.get(c.get("antigen_form_type"), 9),
                len(codes & CONDITIONAL_CODES),
                sum(r.get("severity") == "review" for r in c.get("risks", [])),
                c["candidate_id"])

    ordered = sorted(usable, key=key)
    if not ordered:
        return None, []
    primary, rest = ordered[0], ordered[1:]
    primary_key = key(primary)
    alternates = []
    for c in rest:
        c_key = key(c)
        if c_key[0] > primary_key[0]:
            reason = "完整抗原形式优先于该形式（form priority）"
        elif c_key[1] > primary_key[1]:
            reason = "该备选有更多需特殊设计/核验的问题"
        elif c_key[2] > primary_key[2]:
            reason = "该备选有更多待审阅风险记录"
        else:
            reason = "并列后按确定性次序（candidate_id）排列"
        alternates.append({"candidate_id": c["candidate_id"], "antigen_form_type": c.get("antigen_form_type"),
                           "start": c.get("start"), "end": c.get("end"), "reason_not_primary": reason})
    return primary, alternates


def recommend(protein, candidates, base_flags, entry=None, lab_rules=None):
    """Transparent decision table for one in-scope reference.

    A standard candidate requires positive sequence/topology/boundary evidence
    (a usable candidate already implies sourced exact boundaries); missing
    epitope or contextual-risk literature is recorded as not evaluated and
    never auto-downgrades a candidate with reliable boundaries. Blocking checks
    (coordinate errors, retained native TM/SP/cytoplasmic tail) keep working:
    blocked candidates are excluded before classification.
    """
    entry = entry or {}
    topology = protein.get("topology", "unknown")
    base_codes = {f["code"] for f in base_flags}
    record = {
        "entry_id": entry.get("entry_id", ""), "protein_id": protein.get("protein_id", ""),
        "gene": protein.get("gene", ""), "accession": protein.get("accession", ""),
        "isoform": protein.get("isoform", ""),
        "reference_selection": protein.get("reference_selection"),
        "experimental_isoform_confirmation": (protein.get("reference_selection") or {}).get("experimental_isoform_confirmation", "not_performed"),
        "topology": topology, "location": protein.get("location", "unknown"),
        "scope_status": "in_scope", "scope_reason_codes": [],
        "extension_set": "core",
        "processing_status": "ok", "processing_error": "",
        "screening_recommendation": None, "reason_codes": [], "rationale": "",
        "evidence": {}, "missing_info": [], "not_evaluated": [],
        "candidates_evaluated": len(candidates),
        "primary_candidate_id": "", "alternates": [],
        "lab_rules": {"status": "applied" if lab_rules else "尚未纳入",
                      "note": "未导入实验室经验规则；建议未由实验数据校准" if not lab_rules else ""},
        "functional_status": "not_experimentally_validated",
    }
    disp = entry.get("disposition") or disposition(protein)
    record["scope_status"] = disp.get("scope_status", "in_scope")
    record["scope_reason_codes"] = disp.get("scope_reason_codes", [])
    record["extension_set"] = disp.get("extension_set", "core")
    if record["scope_status"] == "out_of_scope":
        record["rationale"] = "范围外对象：不进入设计评估；这不证明任何工程展示路线不可能"
        record["reason_codes"] = list(record["scope_reason_codes"])
        return record
    lacking = sorted(base_codes & INSUFFICIENT_CODES)
    if lacking:
        record["screening_recommendation"] = "insufficient_evidence"
        record["reason_codes"] = lacking
        record["rationale"] = "身份、拓扑或关键边界等核心信息不足或相互冲突，无法提出可靠建议；这是证据状态，不是生物学失败"
        record["missing_info"] = lacking
        return record
    primary, alternates = pick_primary(candidates)
    if primary is None:
        blocked_codes = sorted({r["code"] for c in candidates for r in c.get("risks", []) if r.get("severity") == "block"} | {c for c in base_codes})
        record["screening_recommendation"] = "no_standard_route"
        record["reason_codes"] = blocked_codes or ["no_supported_continuous_candidate"]
        record["rationale"] = ("在当前定义的候选路线（完整抗原形式 / 有来源的成熟链、GPI 成熟形式、shed 形式或结构域备选）"
                               "与现有证据下没有支持的常规方案；这是路线限制，不表示该蛋白永远不可用")
        record["missing_info"] = ["缺少可用的候选路线支持：", *record["reason_codes"]]
        return record
    codes = sorted({r["code"] for r in primary.get("risks", [])})
    conditional = sorted(set(codes) & CONDITIONAL_CODES)
    if "prediction_requires_annotation_review" in conditional and protein.get("evidence", {}).get("kind") == "curated_annotation":
        # Reviewed entry: prediction-evidence boundaries remain a review-level
        # risk but do not downgrade the class. Unreviewed entries stay conditional.
        conditional.remove("prediction_requires_annotation_review")
    cls = "conditional_candidate" if conditional else "standard_candidate"
    record.update(screening_recommendation=cls, primary_candidate_id=primary["candidate_id"], alternates=alternates)
    record["evidence"] = {
        "boundary_evidence": primary.get("boundary_evidence", {}),
        "topology_basis": {"topology": topology, "location": protein.get("location", "unknown")},
        "reference": {"accession": protein.get("accession"), "isoform": protein.get("isoform"),
                      "evidence": protein.get("evidence", {}), "selection": protein.get("reference_selection")},
        "sequence_sha256": primary.get("reference_sha256", ""),
        "antigen_form_type": primary.get("antigen_form_type"),
        "interval": {"start": primary.get("start"), "end": primary.get("end"), "coordinate_system": "1-based-inclusive"},
        "length": primary.get("length"),
    }
    record["reason_codes"] = conditional if conditional else ["positive_sequence_topology_boundary_evidence"]
    if cls == "standard_candidate":
        record["rationale"] = ("有明确参考序列、适用的拓扑与区间依据，能够提出常规候选；"
                               "这证明候选在计算上成立，不证明表达或激活成功")
    else:
        record["rationale"] = "有具体、有来源的候选，但存在需要特殊设计或重点核验的明确问题：" + ", ".join(conditional)
    # Missing literature/risk records are reported, never converted into failure.
    missing = []
    criteria = {r["kind"]: r for r in primary.get("criteria_review", [])}
    absent = [k for k in CONTEXTUAL_KINDS if criteria.get(k, {}).get("status") == "missing"]
    if absent:
        missing.append("上下文风险记录缺失（未评估，非低风险）: " + ", ".join(absent))
    if any(r["code"] == "epitope_coverage_unknown" for r in primary.get("risks", [])):
        missing.append("无已映射表位：表位保留未评估（完整 ECD 不等于识别保留已被证明）")
    if record["experimental_isoform_confirmation"] != "performed":
        missing.append("实验室实际 isoform 未实验确认；公共 canonical 仅作分析参考，融合装配前必须确认")
    if protein.get("location") != "plasma_membrane" and topology != "secreted":
        missing.append("质膜定位注释未确认")
    record["missing_info"] = missing
    not_evaluated = ["scaffold_assembly（缺真实骨架，不生成融合序列）", "functional_endpoints（四类实验终点均未测定）"]
    if absent:
        not_evaluated.append("contextual_risk_literature（六类上下文风险文献未逐条核读）")
    if primary.get("risk_review_status") == "incomplete":
        not_evaluated.append("risk_review_incomplete")
    record["not_evaluated"] = not_evaluated
    applied = _lab_rules_for(protein, lab_rules)
    if lab_rules is not None:
        record["lab_rules"] = {"status": "applied" if applied else "no_rule_matched", "applied": applied,
                               "note": "" if applied else "已导入实验室经验规则，无匹配条目"}
    for rule in applied:
        code = rule.get("reason_code") or ("lab_rule_" + rule["rule_id"])
        if code not in record["reason_codes"]:
            record["reason_codes"].append(code)
        if rule["effect"] == "downgrade_to_conditional" and cls == "standard_candidate":
            record["screening_recommendation"] = cls = "conditional_candidate"
            record["rationale"] += f"；实验室经验规则 {rule['rule_id']} 调整为条件性候选：{rule['rationale']}"
        elif rule["effect"] == "downgrade_to_no_standard_route" and cls in {"standard_candidate", "conditional_candidate"}:
            record["screening_recommendation"] = "no_standard_route"
            record["rationale"] += f"；实验室经验规则 {rule['rule_id']} 记录该路线在实验室条件下无常规方案：{rule['rationale']}"
    return record


def run_screening(set_dir, outdir, lab_rules_path=None, domain_policy="auto", pilot=False, resume=False):
    """Batch screening bundle: per-item isolation, resume, partial marking."""
    if domain_policy not in {"auto", "full_ecd_only", "full_ecd_and_domains"}:
        raise ValueError("Unknown domain policy")
    definition, proteins = load_set_proteins(set_dir)
    lab_rules = load_lab_rules(lab_rules_path) if lab_rules_path else None
    options = {"domain_policy": domain_policy, "pilot": bool(pilot),
               "lab_rules_sha256": file_hash(lab_rules_path) if lab_rules_path else None}
    code_hash = engine_hash()
    run_key = digest({"set": definition["set_id"], "entries": definition["entries"],
                      "options": options, "engine_sha256": code_hash})
    root = Path(outdir)
    root.mkdir(parents=True, exist_ok=True)
    run_dir = root / run_key[:20]
    attempt = 1
    while run_dir.exists():
        if resume and verify_bundle(run_dir):
            return {"run_dir": str(run_dir.resolve()), "execution": "verified_cache_hit",
                    "summary": read_json(run_dir / "summary.json")}
        if not resume:
            raise FileExistsError("Screening run exists; use --resume or a new output root")
        attempt += 1
        run_dir = root / f"{run_key[:20]}-attempt-{attempt}"
        if verify_bundle(run_dir):
            return {"run_dir": str(run_dir.resolve()), "execution": "verified_cache_hit",
                    "summary": read_json(run_dir / "summary.json")}
    run_dir.mkdir()
    state = {"run_key": run_key, "engine_sha256": code_hash, "version": __version__,
             "started_at": now(), "state": "running", "phases": [], "set_id": definition["set_id"],
             "options": options}

    def event(phase, status, **details):
        item = dict(phase=phase, status=status, timestamp=now(), **details)
        with (run_dir / "events.jsonl").open("a", encoding="utf-8") as handle:
            handle.write(canonical(item) + "\n")
        state["phases"].append(item)
        write_json(run_dir / "state.json", state)

    try:
        write_json(run_dir / "set_definition.json", definition)
        event("set_loaded", "complete", entries=len(definition["entries"]))
        records, all_candidates = [], []
        for entry in definition["entries"]:
            try:
                record, candidates = _screen_entry(entry, proteins, lab_rules, domain_policy)
            except Exception as exc:  # per-item isolation
                record = {"entry_id": entry["entry_id"], "protein_id": entry.get("accession", ""),
                          "gene": entry.get("gene", ""), "accession": entry.get("accession", ""),
                          "isoform": "", "reference_selection": entry.get("reference_selection"),
                          "experimental_isoform_confirmation": "not_performed",
                          "topology": "unknown", "location": "unknown",
                          "scope_status": "unevaluated", "scope_reason_codes": [], "extension_set": "core",
                          "processing_status": "technical_failure", "processing_error": f"{type(exc).__name__}: {exc}",
                          "screening_recommendation": None, "reason_codes": ["technical_failure"],
                          "rationale": "处理过程出现技术失败；这不是生物学结论", "evidence": {},
                          "missing_info": [], "not_evaluated": [], "candidates_evaluated": 0,
                          "primary_candidate_id": "", "alternates": [],
                          "lab_rules": {"status": "尚未纳入" if not lab_rules else "not_applied"},
                          "functional_status": "not_experimentally_validated"}
                candidates = []
            records.append(record)
            all_candidates.extend(candidates)
        event("screening", "complete", records=len(records), candidates=len(all_candidates))
        # Annotate candidate copies with screening context; design records stay intact.
        primary_by_protein = {r["protein_id"]: r for r in records if r.get("primary_candidate_id")}
        for c in all_candidates:
            rec = primary_by_protein.get(c.get("protein_id"))
            c["screening_recommendation"] = rec["screening_recommendation"] if rec else None
            c["is_primary"] = bool(rec and c["candidate_id"] == rec["primary_candidate_id"])
            protein = proteins.get(rec["entry_id"]) if rec else None
            c["junction_notes"] = junction_notes(c, protein or {})
        from .reporting import export_screening_outputs
        summary = export_screening_outputs(run_dir, definition, records, all_candidates,
                                           version=__version__, engine_sha256=code_hash, pilot=pilot)
        event("export", "complete")
        state["state"] = "complete"
        state["set_completeness"] = summary["completeness"]["state"]
        files = {p.name: file_hash(p) for p in sorted(run_dir.iterdir()) if p.is_file() and p.name != "manifest.json"}
        write_json(run_dir / "manifest.json", dict(state, artifacts=files, finished_at=now()))
        return {"run_dir": str(run_dir.resolve()), "execution": "computed", "summary": summary}
    except Exception as exc:
        state["state"] = "failed"
        event("failure", "error", error=str(exc))
        raise


def partition_features(protein):
    """Split features into engine input vs explicitly deferred annotations.

    Only screening uses this; the strict assembly path (harness.run) keeps
    every fuzzy/invalid feature blocking. Deferral never applies to the sole
    candidate-essential annotation of its kind, so genuine coordinate
    problems remain blocking. Deferred annotations are recorded, not dropped.
    """
    from .common import interval, sequence, sourced
    seq = sequence(protein.get("sequence", ""))
    features = protein.get("features", [])
    valid_kinds = set()
    validity = []
    for f in features:
        try:
            interval(f, len(seq))
            if not sourced(f):
                raise ValueError("Feature lacks source/version/evidence kind")
            validity.append(True)
            valid_kinds.add(f.get("kind"))
        except ValueError:
            validity.append(False)
    engine_input, deferred = [], []
    for f, is_valid in zip(features, validity):
        if is_valid:
            engine_input.append(f)
            continue
        kind = f.get("kind", "unknown")
        essential = kind in CANDIDATE_ESSENTIAL_KINDS or (
            kind in {"chain", "gpi_attachment_site"} and protein.get("topology") in CHAIN_ESSENTIAL_TOPOLOGIES)
        if essential and kind not in valid_kinds:
            engine_input.append(f)  # sole essential annotation: must stay blocking
        else:
            deferred.append({"kind": kind, "name": f.get("name", ""), "start": f.get("start"),
                             "end": f.get("end"), "reason": "fuzzy/invalid/unsourced coordinates",
                             "note": "次要/备选注释的坐标问题被显式保留为待审阅，未参与候选生成；原始记录未修改"})
    return engine_input, deferred


def _screen_entry(entry, proteins, lab_rules, domain_policy):
    """Screen one set entry; duplicates reference the first occurrence."""
    if entry.get("duplicate_of"):
        return {"entry_id": entry["entry_id"], "protein_id": entry.get("accession", ""),
                "gene": entry.get("gene", ""), "accession": entry.get("accession", ""), "isoform": "",
                "reference_selection": entry.get("reference_selection"),
                "experimental_isoform_confirmation": "not_performed",
                "topology": "unknown", "location": "unknown",
                "scope_status": "unevaluated", "scope_reason_codes": ["duplicate_input"],
                "extension_set": "core", "processing_status": "ok", "processing_error": "",
                "screening_recommendation": None, "reason_codes": ["duplicate_input"],
                "rationale": "重复输入行：保留记录并指向首次出现，不重复计数",
                "evidence": {}, "missing_info": [], "not_evaluated": [], "candidates_evaluated": 0,
                "primary_candidate_id": "", "alternates": [], "duplicate_of": entry["duplicate_of"],
                "lab_rules": {"status": "尚未纳入" if not lab_rules else "not_applied"},
                "functional_status": "not_experimentally_validated"}, []
    protein = proteins.get(entry["entry_id"])
    if protein is None:
        status = "technical_failure" if entry.get("processing_status") == "failed" else "ok"
        reason = "identity_ambiguous" if entry.get("identity_status") == "ambiguous" else \
            "identity_unresolved" if entry.get("identity_status") != "resolved" else "technical_failure"
        record = {"entry_id": entry["entry_id"], "protein_id": entry.get("accession", "") or entry.get("input", {}).get("value", ""),
                  "gene": entry.get("gene", "") or (entry.get("input", {}) or {}).get("value", ""),
                  "accession": entry.get("accession", ""), "isoform": "",
                  "reference_selection": None, "experimental_isoform_confirmation": "not_performed",
                  "topology": "unknown", "location": "unknown",
                  "scope_status": "unevaluated", "scope_reason_codes": [], "extension_set": "core",
                  "processing_status": status, "processing_error": entry.get("processing_error", ""),
                  "screening_recommendation": None if status == "technical_failure" else "insufficient_evidence",
                  "reason_codes": [reason],
                  "rationale": "身份未解析或存在歧义（已保留全部候选映射，不做猜测）；这不是生物学失败" if status == "ok"
                  else "获取或整理过程失败；这不是生物学结论",
                  "evidence": {"mapping": entry.get("mapping"), "fetch": entry.get("fetch")},
                  "missing_info": [reason], "not_evaluated": [], "candidates_evaluated": 0,
                  "primary_candidate_id": "", "alternates": [],
                  "lab_rules": {"status": "尚未纳入" if not lab_rules else "not_applied"},
                  "functional_status": "not_experimentally_validated"}
        return record, []
    p = dict(protein)
    if domain_policy == "full_ecd_and_domains" or (domain_policy == "auto" and p.get("topology") == "multi_pass"):
        # Multi-pass has no single ECD route; sourced external domains are the
        # only candidate material, so the pipeline explicitly opts in and the
        # resulting domain fragments remain conditional.
        p["candidate_policy"] = "full_ecd_and_domains"
    deferred = []
    if p.get("sequence") and isinstance(p.get("features"), list):
        try:
            p["features"], deferred = partition_features(p)
        except ValueError:
            pass  # invalid sequence: the engine's blocking check handles it
    candidates, base_flags = propose(p, None)
    record = recommend(p, candidates, base_flags, entry=entry, lab_rules=lab_rules)
    if deferred:
        record["deferred_annotations"] = deferred
        record.setdefault("reason_codes", [])
        if "secondary_feature_annotation_deferred" not in record["reason_codes"]:
            record["reason_codes"].append("secondary_feature_annotation_deferred")
        record["missing_info"] = list(record.get("missing_info", [])) + [
            "次要/备选注释坐标待审阅（不影响主候选，装配前必须解决）: "
            + ", ".join(f"{d['kind']}:{d.get('name') or ''}" for d in deferred)]
    return record, candidates
