# v0.5.0: Candidate Decisions, Evidence and Receiver Review

## 目的与主线

面向 **Antibody-Sender → Antigen-Receiver**，回答每个已定义参考“可提出什么片段、为什么、取舍是什么、下一步核验什么”，不是预测 SNIPR 成功率。保留候选引擎、坐标/SP/TM/胞内尾等硬检查；不生成新骨架、不做随机截短、不要求逐蛋白调用 LLM。

1. 固定集合、查询版本、身份和参考序列。公共 canonical 不等于实验 isoform。
2. 规范化适用于该参考的注释；不转移其他 isoform 的坐标。
3. 生成完整连续 ECD、成熟链/GPI 或有来源的域备选。
4. 排除硬阻断候选，比较具体结构与表位取舍，保留主要候选和备选。
5. 分别输出初筛类别、注释支持、路线缺口、接收端待核验项。
6. 用分层人工复核推动规则修订；用真实骨架与四类终点推动实验闭环。

## 不可合并的判断轴

| 轴 | 回答的问题 | 不回答什么 |
|---|---|---|
| scope_status | 是否进入声明的天然表面/分泌扩展设计范围 | 范围外不等于人工工程不可能 |
| screening_recommendation | 常规、条件性、当前无路线、核心信息不足 | 不是功能标签或可信度分数 |
| annotation_support | 参考/边界/拓扑/域的具体来源及 ECO | reviewed 不等于实验验证 |
| route_diagnostic | 无候选是缺资料、路线未支持、阻断还是技术问题 | 不推断永久不能构建 |
| design/assembly_status | 片段一致性、实际骨架与审阅是否满足合同 | 不证明表达、识别或激活 |
| functional endpoints | 表面表达、识别保留、背景、诱导四个实测终点 | 不以一个成功标签合并 |

## 主候选选择合同

只在同一分析参考的非阻断、有序列候选中按以下元组排序，不做权重加总：

1. `annotated_integrity_disruption`：是否有已注释域切割或跨边界二硫键，优先无此已知问题的方案。
2. `mapped_epitope_loss`：是否部分/完全去除已映射表位，优先无此已知丢失的方案。
3. `antigen_form_priority`：完整 ECD → 成熟 GPI → 成熟分泌链 → 有释放依据的 shed 形式 → 域备选。
4. `candidate_id_tiebreak`：纯粹保持可复现，不表达生物学优劣。

若完整形式有结构切割而域备选没有，备选可能成为主要候选，但其省略域/表位和独立折叠问题仍明示；完整形式保留为备选，不悄悄删除。已知表位并不穷尽新抗体的潜在表位，因此不能把只保留一个已知表位作为全面筛抗体的最优设计。该排序仍是待验证的设计政策，不是经过实验校准的模型。

长度、Cys 数、sequon 和上下文风险条数不作排序总分。未知表位是 unknown，不等于没有表位、不等于识别已保留。域切割、已注释二硫键跨边界、已映射表位丢失是条件性问题；无表位注释或未核读风险文献不会自动取消候选。

## 注释证据合同

`annotation_support` 保留每个组件的来源、版本、原始 ECO 记录及映射标签。这里是来源分类而非概率：

- ECO:0000269 → experimental_annotation（该注释标示实验依据，不是本实验室接收端验证）。
- ECO:0000250 → similarity_transfer；ECO:0000305 → curator_inference。
- ECO:0000255/0000259 → sequence_model_inference；ECO:0000256 → automatic_assertion。
- ECO:0000303 → author_statement；ECO:0000313 → imported_information。
- 未收录 ECO → unclassified_eco 并保留原码；缺 ECO 的人工注释 → curated_support_unspecified。
- synthetic_fixture 永远为 synthetic_only；local_experiment 不自动外推原实验上下文。

不进行 ECO 全本体推理，也不把多个来源合成“高置信”。v0.5.0 取消 reviewed/unreviewed 对预测边界的分类特权；设计路线和注释证据独立汇报。因此“standard + prediction-supported boundary”可以存在，仍需核验，不能脱离证据列单独报告“推荐已可靠”。缺失正面序列/拓扑/区间注释仍不能成为常规候选。

依据：[EMBL-EBI UniProt evidence](https://www.ebi.ac.uk/training/online/courses/uniprot-exploring-protein-sequence-and-functional-info/where-does-the-data-come-from/data-evidence/) 说明人工整理可包含实验、计算和转移注释；[UniProt manual](https://prosite.expasy.org/docs/userman.html) 定义逐注释证据及来源记录。两者不为本工具的主候选排序提供生物学验证。

## 路线诊断

- `core_annotation_gap_or_conflict`：身份/序列/拓扑/胞外边界/GPI 成熟边界缺失或冲突，补核心证据。
- `no_supported_route_in_current_annotations`：当前注释与所选策略无可用路线，例如多跨膜只有胞外环而无完整外部域；不是证明不存在自主折叠域。
- `proposed_routes_blocked`：已提议方案被具体硬规则阻断，检查 reason_codes，不能直接绕过。
- `lab_context_route_restriction`：有来源的实验室规则限制当前路线，不外推所有条件。
- `outside_declared_scope`：当前范围外，保留天然定位及单列分母，不自动设计细胞器腔面扩展库。
- `technical_failure` / `duplicate_input` / `supported_candidate_unvalidated` / `not_assessed` 保持独立。

## 接收端与实验闭环

`receiver_review` / TSV 检查的是实际骨架身份、连接/系留方向、识别范围、结构/搭档、加工/连接区问题。明确缺少实际骨架，绝不从 qDis 发送端 PDGFR-LC 或专利默认件拼出 SNIPR。若抗原槽位于 TM 上游，会以片段 C 端系留；II 型天然 N 端系留因此需特别核验，但这里仍不是完成的接头配置。

四终点初始化为 not_measured；按构建、候选、骨架版本、发送端/抗体、批次、重复、刺激条件、门控分母、单位关联。已有 `run` 分支继续负责真实装配、批准和实验合同。批量 `screen` 不猜历史标签，不把表达阳性当诱导阳性。

## 分层核验而非伪造准确率

`manual_review_queue.tsv` 按 scope、拓扑、初筛类别、路线、边界支持组合、结构/表位问题分层，以稳定哈希选每层最多 N 行（默认2）。这是目的性核验清单，不是等概率样本；抽审不能直接计算全体准确率。同一输入版本/参数/代码产生相同队列。空白结果与 pending 审阅不得转为失败或批准。

建议先由独立人员在不依赖软件结论的条件下记录边界来源和构建取舍，再比较差异。此队列含软件分层信息，因此本身不构成盲法实验。规则校准集与独立验证集分离；未来建模按蛋白/同源家族分组，同基因构建不任意跨组。真正准确率需要预先定义统计口径及独立标签，当前不具备。

## 输出与验收

新增三个 TSV：candidate_comparison、receiver_review_plan、manual_review_queue；JSON 内保留完整字段。每个输入保留 route_diagnostic；summary.methodology 提供逐行路线计数和主候选证据/表位计数，分别给分母。新文件加入原运行 manifest 的 SHA256 校验，不改变旧包。

测试合成边界、unknown/mixed ECO、资料缺口、选择取舍、四个未测终点、队列确定性和参数入哈希；另用公开缓存的小批量 pilot 验收输出和原有基线。软件测试不等于生物学准确率。旧全量结果只作历史，不在本版自动刷新。

## 其他 agent 启动与迁移

先读 SKILL、workflow、data-contracts 和本页。首次 `make all-checks-offline`；使用包内脚本而不是 PATH 上可能过时的安装版本。

```bash
python3 scripts/ecd_snipr_cli.py screen --list /private/targets.tsv \
  --cache /private/public-cache --outdir /private/new-pilot --offline --pilot \
  --review-per-stratum 2 --resume
python3 scripts/ecd_snipr_cli.py verify-run --run-dir /private/new-pilot/RUN_HASH
```

无缓存时仅为公共 ID 去掉 `--offline`。生产全量入口是 build-set --query default（不加 --limit）再 screen --set；仅在用户授权全量运行后执行。旧规范化输入若来自 v0.4.1，必须从原始缓存重建后使用，不覆写旧运行。

部署：先备份现有安装，再 `make install-user`；这是一项独立操作，本次代码升级不默认替换安装。旧软件包和原始工作簿保持不动。给 agent 的启动语可用：

> 使用当前 v0.5.0 ecd-snipr-harness，先核实版本与哈希；按 Antibody-Sender → Antigen-Receiver 做指定名单的候选初筛，先 pilot，不自动全量或发布。分别报告初筛类别、逐注释证据、结构/表位取舍、路线缺口和真实骨架待核验项；无骨架不生成融合序列，无实验不填成功标签。交付可对账表格、FASTA、汇总和运行校验。
