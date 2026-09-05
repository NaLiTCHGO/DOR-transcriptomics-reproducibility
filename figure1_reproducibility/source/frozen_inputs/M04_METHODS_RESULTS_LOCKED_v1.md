# M04 跨队列复现与异质性：冻结方法与结果 v1

`MODULE_STATE: PASS_WITH_LIMITATION`

`RUN_ID: 20260814_M04_B1_THREE_CORE_RNASEQ`

## 研究目的与边界

本模块检验三个独立、已完成队列内 QC 和效应估计的人源卵巢颗粒细胞 RNA-seq 队列中，DOR 相对 NOR 的转录效应是否在方向、排序及幅度层面可复现，并量化基因层面异质性。M04 只读取冻结的 gene-level 效应表，不读取或合并样本表达矩阵，不进行跨队列归一化、批次校正或 megamatrix 建模。

## 冻结输入

- GSE274832：6 个独立样本，condition-only `DOR−NOR`。
- GSE193136：12 个独立样本，年龄校正主模型 `~ age_centered + condition` 的 `DOR−NOR`。
- GSE232306：12 个独立样本，因公开 metadata 缺少个体年龄，采用 condition-only `DOR−NOR`。
- 三个输入均来自 GRCh38.p14 + GENCODE v50 raw-to-Salmon/DESeq2 流程。
- 综合键为不带版本号的 Ensembl `gene_id`；基因符号仅用于显示。

输入行数分别为 16,012、16,656 和 16,064。所有输入 `gene_id` 唯一，所有进入综合的 `lfcSE` 为有限正值。三队列共同有效基因集合为 13,993 个。

## 统计方法

### 主分析

1. 在 13,993 个共同基因中，两两计算 DOR−NOR `log2FoldChange` 的 Spearman 相关。
2. 两两计算非零效应方向一致率。
3. 按绝对 Wald z 值选择 top 100、500 和 1,000，计算集合重叠、Jaccard 及重叠基因方向一致率。
4. 队列内候选固定定义为 `padj<0.05 且 |log2FC|>=1`；候选重叠只作描述，不作为外部验证。

### 辅助 meta-analysis

每个共同基因以 `log2FoldChange` 为效应量、`lfcSE` 为标准误，计算逆方差固定效应、Cochran Q、DerSimonian–Laird tau²、随机效应、BH-FDR 和 I²。I² 分为 `<25%`、`25–<50%`、`50–<75%` 和 `>=75%`。由于仅有三个队列，Q、tau² 和 I² 仅作描述。

严格探索性共识规则在运行前冻结为：随机效应 FDR<0.05、|随机效应 log2FC|>=0.5、三个队列全部同向、至少两个队列达到内部候选阈值且 I²<50%。

### 留一队列分析

依次留出每个队列，使用剩余两个队列的逆方差固定效应，与三队列随机效应估计比较 Spearman 排序、方向一致率及效应位移。该 gene-level LOCO 是 M04 的队列影响审计；M06 仍负责主通路结论的正式 leave-one-cohort-out。

## 主要结果

### 全基因效应复现

| 队列对 | Spearman rho | 方向一致率 | 解释 |
|---|---:|---:|---|
| GSE274832 vs GSE193136 | 0.2345 | 57.70% | 弱正相关 |
| GSE274832 vs GSE232306 | 0.2363 | 59.41% | 弱正相关 |
| GSE193136 vs GSE232306 | 0.0491 | 51.71% | 接近无相关 |

top-500 绝对 Wald z 重叠分别为：

- GSE274832 vs GSE193136：54 个，Jaccard 0.0571，重叠基因方向一致率 90.74%。
- GSE274832 vs GSE232306：0 个。
- GSE193136 vs GSE232306：1 个，且该基因方向不一致。

这些结果说明全局 gene-level 可重复性较弱，尤其 GSE193136 与 GSE232306 的效应排序几乎不相关。

### 队列内候选重叠

在共同基因集合中，达到 `FDR<0.05 且 |log2FC|>=1` 的基因数为：

- GSE274832：21 个（该队列全部 30 个候选中有 21 个位于三队列共同集合）。
- GSE193136：19 个（全部 22 个候选中有 19 个位于共同集合）。
- GSE232306：5,717 个（全部 6,381 个候选中有 5,717 个位于共同集合）。

没有基因在三个队列中同时达到该候选阈值。GSE274832 与 GSE232306 重叠 10 个；GSE193136 与 GSE232306 重叠 9 个；GSE274832 与 GSE193136 重叠 0 个。共有 11 个基因在至少两个候选队列中方向相同。共有 4,815/13,993 个基因不论显著性在三个队列中效应方向一致。

GSE232306 的 5,717 个共同集合候选远高于另外两个队列，结合其 phenotype-aligned GC、重复率及 PC1 分离，支持将其视为广泛 cohort-wide signature，而非 5,717 基因的可靠标志物集合。

### Meta-analysis 与异质性

- 随机效应 FDR<0.05：449 个基因。
- I²<25%：4,357 个基因。
- I² 25–<50%：1,306 个基因。
- I² 50–<75%：2,447 个基因。
- I²>=75%：5,883 个基因。

42.0%（5,883/13,993）的共同基因处于 I²>=75% 的高异质性区间，说明大量 gene-level 效应在队列间不稳定。

严格探索性共识规则仅保留两个基因：

| Gene | GSE274832 log2FC | GSE193136 log2FC | GSE232306 log2FC | Random log2FC | Random FDR | I² |
|---|---:|---:|---:|---:|---:|---:|
| FMNL1 | -2.963 | -1.418 | -2.203 | -2.258 | 0.000944 | 35.65% |
| PGAP1 | 1.811 | 1.408 | 2.180 | 1.766 | 0.00000236 | 41.74% |

FMNL1 在 GSE274832 与 GSE232306 达到内部候选阈值；PGAP1 在 GSE193136 与 GSE232306 达到内部候选阈值。二者均为探索性跨队列共识效应，不是经过临床验证的标志物。

### 留一队列稳定性

| 留出队列 | 与全三队列随机效应的 rho | 方向一致率 | 中位绝对效应位移 |
|---|---:|---:|---:|
| GSE274832 | 0.9376 | 93.09% | 0.1552 |
| GSE193136 | 0.9365 | 92.26% | 0.3024 |
| GSE232306 | 0.5127 | 65.61% | 0.3325 |

去除 GSE232306 对综合排序和方向影响最大，证实三队列随机效应结果明显受该队列影响。两个严格探索性共识基因在三种留一情形中方向均保持不变，但候选数过少，不能据此形成生物标志物结论。

## 模块判定

M04 的输入、计算、25 项独立结构审计和 6 张图的视觉 QC 全部通过；不存在样本级跨队列合并。因此技术上满足 DoD。由于三队列上游均有预设限制、GSE232306 的综合影响明显、全局 gene-level 复现较弱且 k=3 异质性估计不精确，最终判定为 `PASS_WITH_LIMITATION`。

该结果不削弱项目主问题，反而支持既定的 phenotype/provenance-aware heterogeneity 设计。下一步应进入 M05，检验 gene-level 异质性是否在预先冻结的通路/模块层面出现收敛。

## 图形及用途

1. `M04_PAIRWISE_EFFECT_SCATTER.png`：展示三组两两全基因效应关系。
2. `M04_REPRODUCIBILITY_MATRICES.png`：汇总 Spearman rho 与方向一致率。
3. `M04_CANDIDATE_OVERLAP.png`：对数尺度展示阈值候选组合，明确不是 biomarker validation。
4. `M04_HETEROGENEITY_DISTRIBUTION.png`：展示 I² 分布及区间计数。
5. `M04_TOP_CROSS_COHORT_EFFECT_HEATMAP.png`：展示 FMNL1 与 PGAP1 的三队列效应。
6. `M04_LOCO_STABILITY.png`：展示留一队列后的综合效应排序和方向稳定性。

## 允许表述

- 三个独立颗粒细胞 RNA-seq 队列在 gene-level 显示弱、且高度异质的 DOR 转录效应复现。
- M04 定量了方向、排序、候选重叠、随机效应和异质性。
- FMNL1 与 PGAP1 可称为“按预设严格规则得到的探索性共识效应”。
- GSE232306 对综合结果影响显著，必须与其未解析混杂限制共同报告。

## 禁止表述

- 不得把 FMNL1、PGAP1 或其余 449 个 meta-FDR 基因称为已验证标志物、诊断基因、因果基因或临床预测因子。
- 不得把 M04 描述为样本级联合分析、独立临床验证或已完成通路收敛。
- 不得省略 GSE232306 潜在混杂、GSE274832 小样本、GSE193136 低映射率/年龄背景和 k=3 异质性不确定性。

## 可追溯性

- 正式脚本：`04_code/python/M04_01_run_cross_cohort_synthesis.py`。
- 独立审计：`04_code/python/M04_02_audit_cross_cohort_results.py`。
- 统一入口：`04_code/powershell/run_m04_repro_heterogeneity.ps1`。
- 输入契约：`02_protocol_and_design/module_contracts/M04_REPRO_HETEROGENEITY.v1.json`。
- 完整结果：本 run 的 `results/`。
- 审计：`results/M04_AUDIT_REPORT.md`、`results/M04_AUDIT_CHECKS.csv`。
