# E-MTAB-391 Provenance Resolution

日期：2026-08-13  
结论：`RESOLVED_WITH_ANALYSIS_LIMITATION`

## 28 与 26 的真实含义

ArrayExpress SDRF 有 56 行，但不是 56 个样本。每个 `Source Name` 出现两行，分别指向 `input.txt` 和 `normalized.txt`，所以共有 28 个唯一 source/hybridization，也与 processed normalized matrix 的 28 列完全一致。

这 28 个 microarray cycle-samples 包括：

- 15 个 DOR cycle-samples；
- 13 个 NOR egg-donor samples；
- 论文明确说明两名 DOR 患者各在两个周期贡献样本；
- 因此总独立患者数为 26：13 DOR + 13 NOR。

旧评估中“发表分析 cohort 是 13 DOR + 13 NOR、仓库多出两条记录”的表述不准确。正确表述是：论文微阵列分析本身使用 15 DOR + 13 NOR 共 28 张芯片，但来自 13 DOR + 13 NOR 共 26 名患者。

## 能否识别哪两个 DOR source 是重复患者

不能。公开 SDRF 只提供 15 个 cycle-level Source Name，没有 patient ID 或 repeated-patient link；论文也只说明“两名 DOR 患者 represented twice”，未公开对应 source 名称。

因此不能安全地把 15 个 DOR 列任意删成 13 个，也不能在不知道配对关系时把全部 28 列当作 28 个独立患者。

## 安全分析角色

E-MTAB-391 可以保留，但不得作为需要精确 patient-level independence 的主效应估计核心队列。推荐角色：

1. `LEGACY_SUPPORT / SENSITIVITY_COHORT`；
2. 描述性 PCA、方向一致性、基因/通路 rank 或 signed-effect sensitivity；
3. 在 Methods/Table 1 明确写 `N_omics=28 cycle-samples; N_independent=26 patients; two unidentified repeated DOR cycles`；
4. 不报告把 28 列当作 28 名独立患者所得的精确显著性为主证据；
5. 主要跨队列结论须在移除 E-MTAB-391 后仍成立。

如果后续必须进行 formal model，可采用预先声明的保守敏感性路线，例如 DOR cycle-level analysis 与 leave-E-MTAB-out analysis 并列，但不能假称重复患者身份已知。

## 结果位置

- `results/E-MTAB-391_SAMPLE_CYCLE_MANIFEST.csv`
- `raw_metadata/E-MTAB-391.sdrf.txt`
- `raw_metadata/E-MTAB-391.idf.txt`
- `raw_metadata/normalized.txt`

