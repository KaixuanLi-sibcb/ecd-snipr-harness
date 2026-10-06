"""Private, offline laboratory evidence intake. No ranking or model training."""

from collections import Counter, defaultdict
from pathlib import Path
from . import __version__
from .common import canonical, digest, file_hash, now, read_json, tsv, write_json
from .contracts import EVIDENCE_FIELDS, evidence_v2
from .harness import engine_hash, verify_bundle
from .ingest import inventory
from .lab_intake import DICTIONARY_VERSION, extract, field_values, nonblank
from .lab_linkage import (assay_records, audited_constructs, fragment_records, identity_conflicts,
                          load_references, proposed_links, reviewed_assay_conflicts)
from .observations import ENDPOINTS

CELL_COLUMNS = "sheet row cell row_id cell_id field raw_header header_row header_cell raw_value raw_formula formula_present formula_attributes shared_formula_anchors cell_type style_index cell_role mapping_basis present_in_source workbook_sha256".split()
FRAGMENT_COLUMNS = "fragment_record_id sheet row source_cell gene accession agid receiver_plasmid receiver_label sequence_status reference_status sequence sequence_digest matches issues full_construct_identity experimental_isoform_confirmation".split()
ASSAY_COLUMNS = "observation_id row_id field record_kind raw_header reported_endpoint raw_value raw_formula formula_present formula_attributes shared_formula_anchors semantic_status construct_id endpoint unit normalized_value eligibility issues source_ref context raw_record_digest".split()


def jsonl(path, records):
    with Path(path).open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(canonical(record) + "\n")


def raw_evidence(cells, timestamp):
    for cell in cells:
        if cell["cell_role"] != "data" or cell["field"] == "unmapped" or not cell["workbook"].lower().endswith(".xlsx"):
            continue
        r = {k: "" for k in EVIDENCE_FIELDS}
        r.update(record_version="2.0", evidence_type="raw_workbook_cell:" + cell["field"],
                 evidence_value=str(cell["raw_value"]) if nonblank(cell["raw_value"]) else "",
                 evidence_status="confirmed_local_xlsx" if nonblank(cell["raw_value"]) else "missing",
                 source_name=cell["workbook"], source_authority_level="local_lab",
                 source_url="local-cell:" + cell["cell_id"], source_version=cell["workbook_sha256"],
                 retrieved_at=timestamp, query=canonical({"sheet": cell["sheet"], "cell": cell["cell"], "header": cell["raw_header"]}),
                 adapter_name="lab_workbook_intake", adapter_version=__version__,
                 raw_response_sha256=cell["workbook_sha256"], confidence="source_cell_presence_only_not_biological_validation",
                 conflict_status="not_checked", failure_mode="" if nonblank(cell["raw_value"]) else "blank_not_negative",
                 license_note="private institutional evidence; local only")
        evidence_v2(r)
        yield r


def run_lab_intake(source, outdir, *, layout_path=None, semantics_path=None, reference_set=None,
                   construct_run=None, resume=False):
    data = inventory(source)
    layout = read_json(layout_path) if layout_path else None
    semantics = read_json(semantics_path) if semantics_path else None
    if semantics and (semantics.get("schema_version") != "1.0" or semantics.get("source_sha256") != data["source_sha256"]):
        raise ValueError("Semantic contract version or workbook hash mismatch")
    extracted = extract(data, layout)
    wanted = {str(v).strip() for r in extracted["rows"] for v in field_values(r, "accession")}
    references, reference_errors, ref_source = load_references(reference_set, wanted)
    constructs, construct_hash = audited_constructs(construct_run)
    inputs = {"source_sha256": data["source_sha256"], "source_name": data["source_name"],
              "layout": layout, "semantics": semantics, "reference_set": ref_source,
              "loaded_reference_digest": digest(references), "reference_errors": reference_errors,
              "construct_manifest_sha256": construct_hash, "engine_sha256": engine_hash(),
              "dictionary_version": DICTIONARY_VERSION}
    run_id = digest(inputs)[:20]
    root = Path(outdir).resolve() / run_id
    if root.exists():
        if resume and verify_bundle(root) and read_json(root / "inputs.json") == inputs:
            return {"run_dir": str(root), "execution": "reused_verified", "summary": read_json(root / "summary.json")}
        raise FileExistsError("Lab bundle exists; use verified --resume or a new output directory")
    root.mkdir(parents=True)
    state = {"state": "running", "run_id": run_id, "kind": "private_lab_evidence",
             "version": __version__, "started_at": now(), "inputs_digest": digest(inputs)}
    write_json(root / "manifest.json", state)
    try:
        fragments = fragment_records(extracted["rows"], references)
        observations, template = assay_records(extracted["rows"], semantics, constructs)
        conflicts = identity_conflicts(fragments) + reviewed_assay_conflicts(observations)
        links, readiness = proposed_links(extracted["rows"], fragments)
        counts_by_sheet = defaultdict(Counter)
        for cell in extracted["cells"]:
            if cell["cell_role"] == "data" and nonblank(cell["raw_value"]):
                counts_by_sheet[cell["sheet"]][cell["field"]] += 1
        routed_readouts = [o for o in observations if o["record_kind"] == "assay_readout" and nonblank(o["raw_value"])]
        readouts = [o for o in routed_readouts if "layout_requires_review" not in o["issues"]]
        eligible = [o for o in observations if o["eligibility"] == "eligible_measurement_not_success_label"]
        definition_groups = defaultdict(list)
        for o in observations:
            if nonblank(o["raw_value"]):
                definition_groups[(o["source_ref"]["sheet"], o["field"], o["raw_header"])].append(o)
        groups = [{"sheet": key[0], "field": key[1], "raw_header": key[2],
                   "nonblank_source_cells": len(values), "layout_quarantined_cells": sum("layout_requires_review" in o["issues"] for o in values),
                   "reported_endpoints": sorted({o["reported_endpoint"] for o in values}),
                   "example_cells": [o["source_ref"]["cell"] for o in values[:3]],
                   "definition_status": "per_cell_review_required_no_batch_inference"}
                  for key, values in definition_groups.items()]
        missingness = Counter(issue for o in readouts for issue in set(o["issues"]))
        summary = {"version": __version__, "source_sha256": data["source_sha256"],
                   "execution_completeness": "partial" if reference_errors else "complete",
                   "interpretation_completeness": "pending_definitions_and_construct_linkage",
                   "counts": dict(extracted["counts"], fragment_records=len(fragments),
                                  valid_fragment_strings=sum(bool(f["sequence"]) for f in fragments),
                                  unique_fragment_strings=len({f["sequence_digest"] for f in fragments if f["sequence_digest"]}),
                                  reference_status=dict(Counter(f["reference_status"] for f in fragments)),
                                  header_routed_readout_cells_nonblank=len(routed_readouts),
                                  layout_quarantined_readout_cells=len(routed_readouts)-len(readouts),
                                  readout_source_cells_nonblank=len(readouts), eligible_measurements=len(eligible),
                                  explicit_self_activation_source_cells=sum(o["field"] == "self_activation" and nonblank(o["raw_value"]) for o in observations),
                                  supplied_audited_construct_ids=len(constructs), potential_conflicts=len(conflicts),
                                  reference_technical_errors=len(reference_errors),
                                  assay_rows_with_link_hints=sum(r["candidate_fragment_records"] > 0 for r in readiness),
                                  link_readiness=dict(Counter(r["link_readiness"] for r in readiness))),
                   "nonblank_cells_by_sheet_and_field": {s: dict(sorted(c.items())) for s, c in counts_by_sheet.items()},
                   "four_endpoints": {e: {"reported_source_cells": sum(o["reported_endpoint"] == e for o in readouts),
                                          "definition_confirmed_source_cells": sum(o["endpoint"] == e and o["semantic_status"].startswith("confirmed") for o in readouts),
                                          "eligible_measurements": sum(o["endpoint"] == e for o in eligible),
                                          "eligible_distinct_constructs": len({o["construct_id"] for o in eligible if o["endpoint"] == e})} for e in ENDPOINTS},
                   "missingness_reason_counts": dict(missingness),
                   "denominators": {"counts": "source cells/rows and fragment strings as individually named; not independent constructs",
                                    "missingness": "nonblank assay-readout cells; reasons overlap, do not sum as unique observations"},
                   "model_readiness": "not_established_no_success_labels_or_training",
                   "limitations": ["source-cell presence is not biological validity", "no joining by gene, AgID or fragment identity alone",
                                   "reference match does not confirm experimental isoform or complete fusion", "no assay pooling or ranking changes",
                                   "unmapped cells retained; layout detection remains provisional", "no new external API calls"]}
        objects = {"inputs": inputs, "raw_inventory": data, "header_sections": extracted["headers"],
                   "row_records": extracted["rows"], "fragments": fragments, "observations": observations,
                   "reference_errors": reference_errors, "conflicts": conflicts, "summary": summary,
                   "proposed_links": links, "row_readiness": readiness,
                   "semantic_review_template": {"schema_version": "1.0", "source_sha256": data["source_sha256"], "reviews": template}}
        for name, value in objects.items():
            write_json(root / (name + ".json"), value)
        tsv(root / "source_cells.tsv", extracted["cells"], CELL_COLUMNS)
        tsv(root / "fragment_reference_audit.tsv", fragments, FRAGMENT_COLUMNS)
        tsv(root / "assay_evidence.tsv", observations, ASSAY_COLUMNS)
        tsv(root / "construct_link_candidates.tsv", links, ["assay_row_id", "fragment_record_id", "fragment_row_id", "matching_fields", "fragment_sequence_digest", "status"])
        tsv(root / "row_readiness.tsv", readiness, ["sheet", "row", "row_id", "gene", "accession", "agid", "receiver_plasmid", "receiver_label",
                                                  "layout_status", "nonblank_fields", "candidate_fragment_records", "candidate_unique_fragment_strings",
                                                  "link_readiness", "actual_receiver_identity", "assay_semantics"])
        tsv(root / "field_coverage.tsv", [{"sheet": s, "field": f, "nonblank_source_cells": n}
                                         for s, fields in counts_by_sheet.items() for f, n in fields.items()],
            ["sheet", "field", "nonblank_source_cells"])
        jsonl(root / "evidence_records_v2.jsonl", raw_evidence(extracted["cells"], state["started_at"]))
        jsonl(root / "source_claims.jsonl", ({"claim_type": "raw_cell", "source_cell": c,
                                             "status": "recorded" if nonblank(c["raw_value"]) else "missing",
                                             "interpretation": "raw_presence_only_not_biological_validity"}
                                            for c in extracted["cells"] if c["cell_role"] == "data"))
        tsv(root / "assay_definition_groups.tsv", groups, ["sheet", "field", "raw_header", "nonblank_source_cells",
                                                         "layout_quarantined_cells", "reported_endpoints", "example_cells", "definition_status"])
        (root / "LAB_EVIDENCE_REPORT.md").write_text(report(summary, reference_errors), encoding="utf-8")
        if file_hash(source) != data["source_sha256"]:
            raise ValueError("Source changed during read-only intake")
        state.update(state="complete", finished_at=now(), artifacts={p.name: file_hash(p) for p in sorted(root.iterdir())
                                                                    if p.is_file() and p.name != "manifest.json"})
        write_json(root / "manifest.json", state)
        if not verify_bundle(root):
            raise ValueError("Produced bundle integrity check failed")
        return {"run_dir": str(root), "execution": "computed", "summary": summary}
    except Exception as exc:
        write_json(root / "manifest.json", dict(state, state="failed", error=str(exc)))
        raise


def report(summary, errors):
    c = summary["counts"]
    lines = ["# 实验数据回收与可解释性审计", "",
             "## 核心结论", "",
             f"已只读处理 {c['source_sheets']} 张表、{c['data_rows']} 个非空数据行；原始单元格全部保留。",
             f"回收到 {c['readout_source_cells_nonblank']} 个非空读数单元格，其中明确标注 self-activation 的 {c['explicit_self_activation_source_cells']} 个。",
             f"另有 {c['layout_quarantined_readout_cells']} 个被上方读数表头覆盖的非空单元格因布局疑点隔离，未计入上述读数。",
             "这些是来源单元格数，跨表可能重复，不是独立构建数或独立实验数。",
             f"符合完整解释合同的定量记录：{c['eligible_measurements']}。未达合同不意味着没有实验或实验失败。",
             f"片段记录 {c['fragment_records']} 条，其中有效片段字符串 {c['valid_fragment_strings']} 条、去重字符串 {c['unique_fragment_strings']} 个；不等于独立完整接收端。",
             "", "## 四类终点分开", "", "|终点|表头明确指向的非空单元格|合同可解释记录|", "|---|---:|---:|"]
    for endpoint, values in summary["four_endpoints"].items():
        lines.append(f"|{endpoint}|{values['reported_source_cells']}|{values['eligible_measurements']}|")
    lines += ["", "BFP/myc/mRuby 联合门控字段尚未自动归为任何功能终点；SNIPR Assay Y/N 和后续血清学记录不转成成败标签。",
              "", "## 片段核对", ""]
    for status, n in c["reference_status"].items():
        lines.append(f"- {status}: {n} 条片段来源记录。")
    lines += ["", "精确匹配只说明片段可定位于所给参考序列；不确认实际 isoform、完整融合、表位保留或信号功能。未匹配可能来自其他 isoform、标签、突变、连接序列或记录错误，不能直接判错。",
              "", "## 当前仍需补齐", "",
              "1. 一个实际完整接收端构建及对应质粒/片段记录；不要用公开专利骨架代替。",
              "2. 读数单位、MFI统计量、门控顺序与分母，以及无刺激/刺激条件。",
              "3. 构建、批次、重复、发送细胞/抗体和时间点的逐记录关系。",
              "", "语义模板 semantic_review_template.json 默认全部 pending；需责任人按真实记录确认，并使用经校验的 construct run，才能产生可解释定量记录。",
              "先用 assay_definition_groups.tsv 按原始表头讨论字段含义，再核对每条记录的构建与实验条件；分组仅方便核对，不自动合并批次。",
              "", "## 审计与限制", "",
              f"技术处理状态：{summary['execution_completeness']}；参考文件错误 {len(errors)} 个。完整读取不代表语义完整。",
              f"潜在标识/读数冲突 {c['potential_conflicts']} 组；保留不同片段，不按同名合并。",
              f"{c['assay_rows_with_link_hints']} 个含读数或状态的来源行找到了片段关联线索。construct_link_candidates.tsv 只供核对，不能据此转移功能结论。",
              "source_cells.tsv 保留重复表头列、未知列、公式缓存与原坐标；raw_inventory.json 保留全部原始单元格及合并区。",
              "本轮未训练模型、未改公共候选排名、未装配猜测的融合序列、未调用外部服务。所有本地输出包含私有数据，不得发布。", ""]
    return "\n".join(lines)
