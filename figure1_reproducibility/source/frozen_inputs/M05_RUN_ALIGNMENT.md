# M05 运行对齐：20260814_M05_B1_MSIGDB2026_1

## 当前状态

`COMPLETE / PASS_WITH_LIMITATION / READY_TO_LOCK`

## 已完成步骤

| 步骤 | 脚本/文件 | 运行结果 | 是否可用 | 输出位置 |
|---|---|---|---|---|
| 00 软件预审 | R/包只读检查 | R 4.5.3；fgsea 1.36.2；msigdbr 26.1.0；clusterProfiler 4.18.4 | 可用；正式算法选 fgsea | 本文 |
| 01 官方资源冻结 | Broad MSigDB v2026.1.Hs Hallmark/Reactome GMT | PASS：50/1,839 行；48,686/898,680 bytes | 可用；运行不再依赖网络 | `03_data/reference/gene_sets/MSigDB_v2026.1_Hs/` |
| 02 映射预审 | M04 meta table + 冻结 GMT | PASS：13,394 唯一 symbol；合格集合 49 Hallmark + 1,016 Reactome | 可用 | `M05_ANALYSIS_PLAN.v1.md` |
| 03 设计冻结 | M05 contract/DoD/FMEA/plan/resource manifest | PASS：rank、算法、多重校正、收敛规则和停止规则已冻结 | 可用 | `02_protocol_and_design/`、M05 根目录 |
| 04 正式脚本保存 | `M05_01_run_ranked_pathway_convergence.R`、`M05_02_audit_pathway_convergence.py`、`run_m05_pathway_convergence.ps1` | 已保存并登记 40–42；尚待语法检查 | 待运行 | `04_code/` |
| 05 首次语法检查 | R `parse()` | FAIL：仅资源审计输出中的相对路径正则存在 `\\[` 转义错误；尚未启动 fgsea | 结果未产生；科学设计不受影响 | 已改为固定前缀截取，待复查 |
| 06 语法复查 | R/Python/PowerShell parser | PASS：三类脚本全部通过 | 可运行 | `04_code/` |
| 07 首次实际运行 | 统一 PowerShell 入口 | FAIL：rank/resource audit 已完成，但 fgsea 默认 Windows 子进程加载用户 `.Rprofile` 的 `.Last/savehistory`，退出时报错；未产生 NES 结果 | 仅预审文件可用，未锁定 | 已固定 `SerialParam`/`nproc=1`，不修改用户环境 |
| 08 第二次实际运行 | 串行 fgsea | FAIL：第一个 fgsea 计算后，结果列契约使用旧接口 `nMoreExtreme`，当前 fgsea 1.36.2 实际返回 `log2err`；未写出 NES 表 | 未产生可接受结果 | 已按当前锁定版本改为 `log2err` 并只排序实际存在列 |
| 09 第三次实际运行 | 串行 fgsea + pathway synthesis | FAIL：全部 3×2 fgsea 已计算并写出临时总表，但 data.table 在综合阶段未解析 `..collection` 父作用域语法 | 临时 fgsea 表不锁定 | 已改用与列名不冲突的 `collection_name/pathway_name/cohort_name`；正式整批重跑 |
| 10 第四次实际运行 | 修订后的统一入口 | PASS_WITH_LIMITATION：R 计算退出码 0，Python 审计 26/26 PASS | 数值结果可用，待视觉 QC | 本 run `results/` |
| 11 首次六图视觉 QC | 人工检查 | 5 图 PASS；Reactome top-30 长标签挤压热图且标题裁切，FAIL | 数值可用；该图不可锁定 | 已增加标签换行、15-inch 宽度及短标题 |
| 12 fgsea 数值警告审计增强 | 正式 R/Python 脚本 | 记录 eps=1e-10 下限和 `log2err` 缺失计数；禁止声称更小精确 p 值 | 待重跑确认 | 新增 `M05_FGSEA_NUMERICAL_WARNING_AUDIT.csv` |
| 13 第五次实际运行 | 增强后的正式入口 | FAIL：数值与警告审计表已生成，但 Reactome 标签换行把 30 个标签展开为 36 行，ggplot 拒绝作图 | 数值表不锁定，审计未运行 | 已改为逐个 pathway label 单独 `strwrap` |
| 14 最终正式运行 | 完整修订后的统一入口 | PASS_WITH_LIMITATION：退出码 0，耗时 0.35 分钟 | 可用 | 本 run `results/`、`RUN_COMPLETE.txt` |
| 15 独立结构/规则审计 | `M05_02_audit_pathway_convergence.py` | PASS_WITH_LIMITATION；28/28 检查通过 | 可用 | `results/M05_AUDIT_REPORT.md` |
| 16 最终六图视觉 QC | 人工检查 | PASS；Reactome 热图标题、长标签、数值和色标均可读 | 可用 | `VISUAL_QC_REVIEW.PASS.txt`、`results/plots/` |
| 17 方法结果冻结文本 | `M05_METHODS_RESULTS_LOCKED_v1.md` | PASS：方法、数值、资源、失败修复、限制和允许/禁止表述已写清 | 可用 | 本 run 根目录 |

## 下一步

复制完整轻量结果到 `06_locked_results/modules/M05_PATHWAY_CONVERGENCE/v1/`，更新模块登记表和总进度，然后开放 M06 正式 leave-one-cohort-out。
