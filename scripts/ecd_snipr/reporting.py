"""Screening summary outputs: tables, summary.json, PI summary and figure.

All ratios carry explicit denominators. The figure is dependency-free SVG
labelled 候选设计覆盖 / 初筛建议 (never experimental success rates); pilot runs
are labelled pilot. XLSX display may be added downstream, but TSV/JSON remain
the computational formats.
"""

from collections import Counter
from pathlib import Path
from .common import now, tsv, write_json
from .harness import CANDIDATE_COLUMNS
from .screening import SCREENING_COLUMNS, CLASSES

CANDIDATE_TSV_COLUMNS = CANDIDATE_COLUMNS + ["screening_recommendation", "is_primary", "junction_notes"]

TOPOLOGY_ORDER = ["type_i", "type_ii", "gpi", "multi_pass", "secreted", "intracellular", "unknown"]

CLASS_LABELS = {
    "standard_candidate": "常规候选",
    "conditional_candidate": "条件性候选",
    "no_standard_route": "无常规路线",
    "insufficient_evidence": "证据不足",
}


def _rows_for_screening(records, entries, candidates):
    by_entry = {e["entry_id"]: e for e in entries}
    primary = {c["candidate_id"]: c for c in candidates}
    rows = []
    for r in records:
        entry = by_entry.get(r["entry_id"], {})
        disp = entry.get("disposition") or {}
        p = primary.get(r.get("primary_candidate_id") or "", {})
        risks = sorted({x["code"] for x in p.get("risks", [])}) if p else []
        rows.append({
            "entry_id": r["entry_id"], "protein_id": r.get("protein_id", ""),
            "gene": r.get("gene", ""), "accession": r.get("accession", ""), "isoform": r.get("isoform", ""),
            "reference_basis": (r.get("reference_selection") or {}).get("basis", ""),
            "experimental_isoform_confirmation": r.get("experimental_isoform_confirmation", ""),
            "input_row": (entry.get("input") or {}).get("row", ""),
            "input_value": (entry.get("input") or {}).get("value", ""),
            "duplicate_of": r.get("duplicate_of", ""),
            "membrane_set_membership": disp.get("membrane_set_membership", ""),
            "natural_cell_surface_target": disp.get("natural_cell_surface_target", ""),
            "scope_status": r.get("scope_status", ""), "scope_reason_codes": ",".join(r.get("scope_reason_codes", [])),
            "extension_set": r.get("extension_set", ""),
            "topology": r.get("topology", ""), "location": r.get("location", ""),
            "screening_recommendation": r.get("screening_recommendation") or "",
            "primary_candidate_id": r.get("primary_candidate_id", ""),
            "primary_antigen_form": p.get("antigen_form_type", ""),
            "primary_start": p.get("start", ""), "primary_end": p.get("end", ""), "primary_length": p.get("length", ""),
            "reason_codes": ",".join(r.get("reason_codes", [])),
            "rationale": r.get("rationale", ""),
            "main_risks": ",".join(risks),
            "missing_info": " | ".join(r.get("missing_info", [])),
            "not_evaluated": " | ".join(r.get("not_evaluated", [])),
            "deferred_annotations": ";".join(f"{d['kind']}:{d.get('name') or ''}[{d.get('start')}-{d.get('end')}]"
                                             for d in r.get("deferred_annotations", [])),
            "lab_rules_status": (r.get("lab_rules") or {}).get("status", ""),
            "processing_status": r.get("processing_status", ""),
            "processing_error": r.get("processing_error", ""),
        })
    return rows


def build_summary(definition, records, candidates, version, engine_sha256, pilot):
    """Counts never mix units: genes, references/isoforms and candidates are separate."""
    evaluated = [r for r in records if r.get("screening_recommendation") in CLASSES]
    in_scope = [r for r in evaluated if r.get("scope_status") == "in_scope"]
    core = [r for r in in_scope if r.get("extension_set") == "core"]
    classes = Counter(r["screening_recommendation"] for r in in_scope)
    classes_core = Counter(r["screening_recommendation"] for r in core)
    out_of_scope = [r for r in records if r.get("scope_status") == "out_of_scope"]
    failures = [r for r in records if r.get("processing_status") == "technical_failure"]
    unresolved = [r for r in records if r.get("screening_recommendation") == "insufficient_evidence"
                  and any(c in {"identity_unresolved", "identity_ambiguous"} for c in r.get("reason_codes", []))]
    duplicates = [r for r in records if r.get("duplicate_of")]
    genes = {r["gene"] for r in records if r.get("gene")}
    references = {(r["accession"], r["isoform"]) for r in records if r.get("accession") and r.get("duplicate_of") is None and r["accession"]}
    denom = len(in_scope)

    def ratio(n, d, definition_text):
        return {"numerator": n, "denominator": d, "value": (n / d) if d else None,
                "denominator_definition": definition_text}

    summary = {
        "report_type": "候选设计覆盖 / 初筛建议（computational screening coverage, not experimental success rates）",
        "pilot": bool(pilot),
        "generated_at": now(), "version": version, "engine_sha256": engine_sha256,
        "set": {"set_id": definition["set_id"], "database": definition["database"]["name"],
                "release": definition["database"]["release"], "query": definition.get("query"),
                "inclusion_rules": definition.get("inclusion_rules", []),
                "created_at": definition.get("created_at")},
        "completeness": {"state": definition.get("completeness", "partial"),
                         "note": "partial 表示重跑可能增加数据（获取失败/截断）；绝不是全量完成"},
        "counts": {
            "input_rows": definition["counts"]["input_rows"],
            "duplicate_input_rows": len(duplicates),
            "resolved_identities": definition["counts"]["resolved"],
            "unresolved_identities": definition["counts"]["unresolved"],
            "genes_unique": len(genes),
            "references_unique": len(references),
            "candidates_total": len(candidates),
            "candidates_usable": sum(1 for c in candidates if c.get("sequence") and c.get("design_status") != "blocked"),
            "evaluated_references": denom,
            "core_references": len(core),
            "secreted_extension_references": denom - len(core),
            "screening_classes_in_scope": {k: classes.get(k, 0) for k in CLASSES},
            "screening_classes_core_only": {k: classes_core.get(k, 0) for k in CLASSES},
            "out_of_scope": len(out_of_scope),
            "out_of_scope_reasons": dict(Counter(c for r in out_of_scope for c in r.get("scope_reason_codes", []))),
            "technical_failures": len(failures),
            "identity_unresolved_or_ambiguous": len(unresolved),
        },
        "ratios": {
            "standard_candidate_share": ratio(classes.get("standard_candidate", 0), denom,
                                              "纳入设计评估的唯一参考数（核心+分泌扩展；不含范围外/技术失败/重复行）"),
            "conditional_candidate_share": ratio(classes.get("conditional_candidate", 0), denom,
                                                 "同上"),
            "any_candidate_share": ratio(classes.get("standard_candidate", 0) + classes.get("conditional_candidate", 0),
                                         denom, "同上"),
            "no_standard_route_share": ratio(classes.get("no_standard_route", 0), denom, "同上"),
            "insufficient_evidence_share": ratio(classes.get("insufficient_evidence", 0), denom, "同上"),
            "secreted_extension_share": ratio(denom - len(core), denom,
                                              "分泌扩展对象占纳入评估参考的比例，不计入核心分母"),
            "identity_resolution_share": ratio(definition["counts"]["resolved"], definition["counts"]["input_rows"],
                                               "输入行（含重复行）中解析到唯一 accession 的比例"),
        },
        "topology_by_class": {t: {k: 0 for k in CLASSES} for t in TOPOLOGY_ORDER},
        "interpretation": ("这些是确定性的、有依据的计算初筛建议与候选设计覆盖统计，不是实验成功率。"
                           "所有候选 functional_status=not_experimentally_validated。"),
    }
    for r in in_scope:
        t = r.get("topology", "unknown")
        if t not in summary["topology_by_class"]:
            summary["topology_by_class"][t] = {k: 0 for k in CLASSES}
        summary["topology_by_class"][t][r["screening_recommendation"]] += 1
    return summary


def pi_summary(summary, records):
    """Concise PI-facing summary; leads with research conclusions, not counts of approvals."""
    c = summary["counts"]
    r = summary["ratios"]
    classes = c["screening_classes_in_scope"]
    special = Counter()
    for rec in records:
        if rec.get("screening_recommendation") == "conditional_candidate":
            for code in rec.get("reason_codes", []):
                special[code] += 1
    no_route = Counter()
    for rec in records:
        if rec.get("screening_recommendation") == "no_standard_route":
            for code in rec.get("reason_codes", []):
                no_route[code] += 1
    lines = [
        "# 初筛建议汇总（PI_SUMMARY）", "",
        f"**结论先行**：在 {c['evaluated_references']} 个纳入设计评估的参考中，"
        f"{classes.get('standard_candidate', 0)} 个有常规候选方案，{classes.get('conditional_candidate', 0)} 个需要特殊设计或重点核验，"
        f"{classes.get('no_standard_route', 0)} 个在当前候选路线下没有常规方案，{classes.get('insufficient_evidence', 0)} 个核心资料不足。"
        f"另有 {c['out_of_scope']} 个范围外对象、{c['technical_failures']} 个技术失败。"
        "以上为有依据的计算初筛建议（候选设计覆盖），不是实验成功率；所有候选均未实验验证。", "",
        f"> 运行性质：{'PILOT（小批量验收，不代表全量结果）' if summary['pilot'] else '正式批量运行'}；"
        f"集合完整性：{summary['completeness']['state']}（{summary['completeness']['note']}）", "",
        "## 处理了什么", "",
        f"- 集合：{summary['set']['set_id']}；数据库：{summary['set']['database']}（版本：{summary['set']['release']}）",
        f"- 纳入规则：" + "；".join(summary["set"]["inclusion_rules"]),
        f"- 输入行 {c['input_rows']}（重复行 {c['duplicate_input_rows']}，保留不丢行）；解析到唯一 accession：{c['resolved_identities']}；"
        f"唯一 gene：{c['genes_unique']}；唯一参考/isoform：{c['references_unique']}；候选片段：{c['candidates_total']}（可用 {c['candidates_usable']}）",
        f"- 分析参考选择：数据库 canonical/显式 accession 作为分析参考（依据逐条记录）；实验室实际 isoform 均未实验确认，装配前必须确认。",
        "",
        "## 多少对象能提出候选", "",
        f"- 常规候选 standard_candidate：{classes.get('standard_candidate', 0)}"
        f"（{r['standard_candidate_share']['numerator']}/{r['standard_candidate_share']['denominator']}，分母：{r['standard_candidate_share']['denominator_definition']}）",
        f"- 条件性候选 conditional_candidate：{classes.get('conditional_candidate', 0)}（同上分母）",
        f"- 合计能提出候选：{r['any_candidate_share']['numerator']}/{r['any_candidate_share']['denominator']}",
        f"- 分泌扩展对象 {c['secreted_extension_references']} 个单独统计，不混入核心分母（核心参考 {c['core_references']} 个）。",
        "",
        "## 特殊设计（条件性候选）的主要原因", "",
    ]
    if special:
        lines += [f"- {code}: {n} 个对象" for code, n in special.most_common()]
    else:
        lines.append("- 无条件性候选。")
    lines += ["", "## 哪些在当前路线下没有常规方案（路线限制，非“永远不可用”）", ""]
    if no_route:
        lines += [f"- {code}: {n} 个对象" for code, n in no_route.most_common()]
    else:
        lines.append("- 无。")
    lines += [
        "", "## 目前覆盖不到的范围", "",
        f"- 范围外对象 {c['out_of_scope']} 个：" +
        ("；".join(f"{k} {v} 个" for k, v in c["out_of_scope_reasons"].items()) or "无"),
        "- 细胞器膜/管腔面不当作细胞表面；分泌蛋白仅作扩展集；多跨膜不拼接胞外环，仅有来源的外部域可作条件性方案。",
        "",
        "## 多少是资料不足", "",
        f"- 核心资料不足 insufficient_evidence：{c['screening_classes_in_scope'].get('insufficient_evidence', 0)} 个"
        f"（另 {c['identity_unresolved_or_ambiguous']} 个身份未解析/歧义，输入行保留）。",
        "- 缺失的是身份/拓扑/边界证据，不是生物学失败结论。",
        "",
        "## 尚未实验验证的内容", "",
        "- 所有候选：表达、识别保留、背景激活、诱导响应四类终点均未测定（functional_status=not_experimentally_validated）。",
        "- 实验室实际 isoform 未确认；真实 SNIPR 骨架未提供，本轮不生成任何融合序列，仅交付候选片段 FASTA 与接入说明。",
        "- 缺失的表位/上下文风险文献不自动降级已有可靠边界的候选，相应条目标为“未评估”。",
        "",
        f"生成时间：{summary['generated_at']}；版本：{summary['version']}；引擎 SHA256：{summary['engine_sha256'][:16]}…",
    ]
    return "\n".join(lines) + "\n"


def _bar(x, y, width, height, color, label_left, label_right, max_width):
    return (f'<text x="10" y="{y + height - 5}" font-size="12" fill="#1a1a1a">{label_left}</text>'
            f'<rect x="220" y="{y}" width="{width}" height="{height}" fill="{color}" rx="2"/>'
            f'<text x="{225 + width}" y="{y + height - 5}" font-size="12" fill="#1a1a1a">{label_right}</text>')


def figure_svg(summary):
    """Dependency-free SVG: set disposition and topology x class. Source data in JSON."""
    c = summary["counts"]
    classes = c["screening_classes_in_scope"]
    colors = {"standard_candidate": "#2e7d32", "conditional_candidate": "#f9a825",
              "no_standard_route": "#c62828", "insufficient_evidence": "#6a6a6a"}
    panel_a = [
        ("输入条目（含重复/未解析行）", c["input_rows"], "#455a64"),
        ("身份解析到唯一 accession", c["resolved_identities"], "#455a64"),
        ("纳入设计评估（核心+扩展）", c["evaluated_references"], "#1565c0"),
        ("常规候选 standard_candidate", classes.get("standard_candidate", 0), colors["standard_candidate"]),
        ("条件性候选 conditional_candidate", classes.get("conditional_candidate", 0), colors["conditional_candidate"]),
        ("无常规路线 no_standard_route", classes.get("no_standard_route", 0), colors["no_standard_route"]),
        ("证据不足 insufficient_evidence", classes.get("insufficient_evidence", 0), colors["insufficient_evidence"]),
        ("范围外 out_of_scope", c["out_of_scope"], "#8e24aa"),
        ("技术失败 technical_failure", c["technical_failures"], "#b0b0b0"),
    ]
    max_a = max([n for _, n, _ in panel_a] + [1])
    parts = []
    y = 56
    for label, n, color in panel_a:
        width = max(2, int(420 * n / max_a)) if n else 0
        parts.append(_bar(0, y, width, 18, color, label, str(n), 420))
        y += 26
    panel_a_end = y + 10
    # Panel B: topology x class stacked bars
    tops = [(t, v) for t, v in summary["topology_by_class"].items() if sum(v.values())]
    max_b = max([sum(v.values()) for _, v in tops] + [1])
    yb = panel_a_end + 40
    parts.append(f'<text x="10" y="{yb - 12}" font-size="13" font-weight="bold" fill="#1a1a1a">拓扑 × 初筛类别（纳入评估对象）</text>')
    for topology, values in tops:
        x = 220
        total = sum(values.values())
        parts.append(f'<text x="10" y="{yb + 14}" font-size="12" fill="#1a1a1a">{topology}</text>')
        for cls in CLASSES:
            n = values.get(cls, 0)
            if not n:
                continue
            width = max(2, int(420 * n / max_b))
            parts.append(f'<rect x="{x}" y="{yb}" width="{width}" height="18" fill="{colors[cls]}"/>')
            x += width
        parts.append(f'<text x="{x + 5}" y="{yb + 14}" font-size="12" fill="#1a1a1a">{total}</text>')
        yb += 26
    legend_y = yb + 16
    legend = "".join(
        f'<rect x="{20 + (i % 2) * 330}" y="{legend_y + (i // 2) * 20}" width="12" height="12" fill="{colors[k]}"/>'
        f'<text x="{38 + (i % 2) * 330}" y="{legend_y + (i // 2) * 20 + 11}" font-size="11" fill="#1a1a1a">{k}（{CLASS_LABELS[k]}）</text>'
        for i, k in enumerate(CLASSES))
    title = "候选设计覆盖（初筛建议）" + (" · PILOT 小批量验收" if summary["pilot"] else "")
    height = legend_y + 62
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="900" height="{height}" font-family="Helvetica, Arial, sans-serif">'
        f'<rect width="900" height="{height}" fill="#ffffff"/>'
        f'<text x="10" y="24" font-size="16" font-weight="bold" fill="#1a1a1a">{title}</text>'
        f'<text x="10" y="42" font-size="11" fill="#555">确定性计算初筛建议，不是实验成功率；集合完整性：{summary["completeness"]["state"]}；'
        f'生成时间 {summary["generated_at"]}</text>'
        + "".join(parts) + legend +
        (f'<text x="640" y="24" font-size="20" font-weight="bold" fill="#c62828">PILOT</text>' if summary["pilot"] else "")
        + "</svg>")


def figure_source_data(summary):
    c = summary["counts"]
    return {
        "figure": "候选设计覆盖 / 初筛建议",
        "pilot": summary["pilot"],
        "not": "实验成功率",
        "set_disposition": {
            "input_rows": c["input_rows"], "resolved_identities": c["resolved_identities"],
            "evaluated_references": c["evaluated_references"],
            "screening_classes_in_scope": c["screening_classes_in_scope"],
            "out_of_scope": c["out_of_scope"], "technical_failures": c["technical_failures"],
        },
        "topology_by_class": summary["topology_by_class"],
        "denominator_definitions": {k: v["denominator_definition"] for k, v in summary["ratios"].items()},
    }


def export_screening_outputs(run_dir, definition, records, candidates, version, engine_sha256, pilot):
    run_dir = Path(run_dir)
    summary = build_summary(definition, records, candidates, version, engine_sha256, pilot)
    write_json(run_dir / "screening.json", records)
    write_json(run_dir / "candidates.json", candidates)
    write_json(run_dir / "summary.json", summary)
    write_json(run_dir / "screening_overview_data.json", figure_source_data(summary))
    tsv(run_dir / "protein_screening.tsv", _rows_for_screening(records, definition["entries"], candidates), SCREENING_COLUMNS)
    tsv(run_dir / "candidate_plan.tsv", candidates, CANDIDATE_TSV_COLUMNS)
    with (run_dir / "candidate_fragments.fasta").open("w", encoding="utf-8") as handle:
        for c in candidates:
            if c.get("sequence") and c.get("design_status") != "blocked":
                role = "PRIMARY" if c.get("is_primary") else "ALTERNATE"
                handle.write(f">{c['candidate_id']} {c['protein_id']}:{c['start']}-{c['end']} "
                             f"form={c['antigen_form_type']} {role} UNVALIDATED\n{c['sequence']}\n")
    (run_dir / "screening_overview.svg").write_text(figure_svg(summary), encoding="utf-8")
    (run_dir / "PI_SUMMARY.md").write_text(pi_summary(summary, records), encoding="utf-8")
    return summary
