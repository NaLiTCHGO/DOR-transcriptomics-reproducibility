# DOR 纯生信项目最终总览

`FINAL_STATE: GO_FULL / M00-M06 COMPLETE / M09 REPORTING FREEZE`

生成日期：2026-08-15（Asia/Shanghai）

## 一、最终结论

- 唯一正式决策：`GO_FULL`。
- Project Value Score：85.1/100。
- Evidence Readiness Score：86.5/100（项目展示四舍五入为 86.5）。
- 本项目完成了三个核心人源卵巢颗粒细胞 RNA-seq 队列的 provenance、QC、队列内效应、跨队列异质性、通路收敛和三轮通路 LOCO。

## 二、冻结队列

| Accession | N_recruited/context | N_omics | N_independent | DOR/NOR | 角色 |
|---|---|---:|---:|---:|---|
| GSE274832 | 147 clinical participants | 6 | 6 | 3/3 | Core RNA-seq；metadata template conflict disclosed |
| GSE193136 | 24 recruited | 12 | 12 | 6/6 | Core RNA-seq；低 mapping、年龄校正 |
| GSE232306 | 60 recruited | 12 | 12 | 6/6 | Core RNA-seq；潜在技术/临床混杂，screening-only |
| E-MTAB-391 | 28 cycle-samples | 28 | 26 | 14/14 | Legacy sensitivity；非主推断队列 |

## 三、主要结果

1. M04：13,993 个共同基因；两两效应 rho 0.0491–0.2363；随机效应 FDR<0.05 为 449；I²>=75% 为 5,883。FMNL1/PGAP1 仅为探索性共识效应。
2. M05：使用 Human MSigDB v2026.1.Hs；49 Hallmark、1,016 Reactome 进入正式 fgsea；严格收敛为 8 Hallmark + 47 Reactome。GSE274832/GSE193136 的 pathway ranking 优于 gene-level，但含 GSE232306 的比较并未全面改善。
3. M06：三轮均只用两个保留队列选择、第三队列留出验证；无阈值重调。Universal core 为 1 Hallmark + 7 Reactome，主紧凑结果为负向富集的 `HALLMARK_P53_PATHWAY`。
4. 留出 GSE232306 时 Hallmark 严格复现仅 1/14，因此结论必须写成“窄的内部条件稳健性”，不能泛化为全局通路一致。

## 四、允许与禁止表述

允许：弱 gene-level 复现、高异质性、队列依赖的 pathway convergence、窄的三队列内部 LOCO core。

禁止：外部/临床/前瞻性验证；诊断或治疗性能；负 NES 等于通路功能性抑制；因果机制；把候选基因或 8 条通路称为已验证生物标志物。

## 五、可追溯路径

- 日常恢复入口：`PROJECT_PROGRESS.md`。
- 正式代码：`04_code/`；M06 为 `python/M06_01_run_pathway_loco.py`、`python/M06_02_audit_pathway_loco.py`、`powershell/run_m06_pathway_loco.ps1`。
- 锁定结果：`06_locked_results/modules/`。
- M04/M05/M06 冻结方法结果分别位于相应 `v1/*_METHODS_RESULTS_LOCKED_v1.md`。
- 最终 Excel：`DOR_Project_Final_Overview_20260815.xlsx`。
- Claim 对照：`M09_CLAIM_TRACEABILITY.csv`。

## 六、剩余非阻塞事项

- 提交前再次进行最新文献/预印本 collision refresh。
- 受控数据 PRJCA007454 未获授权，不在本项目主结果中。
- GSE232306 个体年龄/完整批次不可得，限制不能通过当前公共数据消除。
- 没有大型同组织独立临床队列，因此不主张外部验证。
