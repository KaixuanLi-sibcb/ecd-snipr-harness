# 公开参考架构与实际实验骨架

核对日期：2026-09-16。机器可读记录：[public_reference_architecture.json](public_reference_architecture.json)。该记录是**参考资料，不是可装配的骨架配置**。

## 已核对的文献描述，不等于性能证明

- Addgene #79127 对应 `pHR_pGK_LaG17_synNotch_Gal4VP64`，识别端是抗 GFP 的 LaG17 纳米抗体；它不是插入了任意靶抗原 ECD 的 Ag-SNIPR。LaG17 也不要与 LAG3 抗原混淆。[Addgene 官方记录](https://www.addgene.org/79127/)
- 专利描述接收侧 myc、mouse Notch1 来源模块、Gal4-VP64；主实施例使用 K562，成像另有 HEK293T 接收细胞实例，不能把一个实施例的宿主推广成所有配置。[CN114437232B](https://patents.google.com/patent/CN114437232B/zh)
- BFP 与 hCD8 是不同报告上下文，不据此认定它们都在 #79127 同一载体上。[WO2022095916A1](https://patents.google.com/patent/WO2022095916A1/en)
- SEQ ID NO.5/6 对应发送侧 PDGFR–LC 模块及编码序列，不能当成接收侧 SNIPR 核心。[CN114437232B](https://patents.google.com/patent/CN114437232B/zh)

## 在本项目中怎么用

公开参考回答“可能参考什么架构”；实际构建资料回答“我们到底使用什么”。两者分别建档，当前一致性标记为 `unconfirmed`。

当前研究配置仍是 **Antibody-Sender -> Antigen-Receiver**。上述公开资料不直接证明角色互换后的表达、背景或响应，也不能自动确定实际插入位点、linker、标签位置和完整序列。

因此，不因实际骨架未确认而停止 ECD 注释、完整结构域候选、表位范围和风险分析。仅把依赖实际连接配置的结论保持待确认，阻止最终融合序列导出。

本轮核对的是公开页面及专利相关段落，没有下载、逐碱基复核 Addgene 全序列或与实验室载体做序列比对。后续应先明确参考版本，再与一份实际完整构建对应，而不是从发送端 SEQ ID 拼接一个假定接收端。

## 再审查后的限定

CN 与 WO 是同一优先权来源的专利家族，不能计算为两个独立实验重复。专利、质粒目录、原始论文、实验室当前构建四类来源分别记录，不把其中一类当作另一类。

也不能反向断言“专利接收端只能使用抗体”：其质粒描述还包括 ACE2 类识别模块。这仍不等于验证了本项目的批量 Antibody-Sender -> Antigen-Receiver，也不保证任意抗原 ECD 的适配性。依据位于 CN 文本“定量展示系统所需要的质粒”及实施例 1，而不是由发送侧 LC 模块推导。

受体文献中的识别模块与调控性胞外模块是不同概念；synNotch 与 SNIPR 不能仅凭名称互换。该参考目录的核对范围是公开描述与来源关系，性能外推状态为 `not_established`。
