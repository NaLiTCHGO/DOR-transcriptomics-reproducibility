# M06 运行对齐：20260814_M06_B1_PATHWAY_LOCO

## 当前状态

`COMPLETE / PASS_WITH_LIMITATION / LOCK_READY`

## 已完成步骤

| 步骤 | 脚本/文件 | 运行结果 | 是否可用 | 输出位置 |
|---|---|---|---|---|
| 00 M05 接口确认 | 只读检查 locked M05 convergence/pairwise/strict/leading-edge 表 | PASS：49 Hallmark、1,016 Reactome、55 frozen strict pathways | 可用 | `06_locked_results/modules/M05_PATHWAY_CONVERGENCE/v1/` |
| 01 M06 预审 | 临时只读规则核算 | 三次留出均有严格 pathway 复现；仅用于确认设计可计算，不锁定数值 | 不作为正式结果 | 本文 |
| 02 设计冻结 | M06 contract/DoD/FMEA/plan | PASS：pair selection、held-out replication、universal core 和停止规则已冻结 | 可用 | `02_protocol_and_design/`、M06 根目录 |
| 03 正式脚本保存 | `M06_01_run_pathway_loco.py`、`M06_02_audit_pathway_loco.py`、`run_m06_pathway_loco.ps1` | 已保存并登记 43–45 | 可复现 | `04_code/` |
| 04 语法检查 | `python -m py_compile`；PowerShell AST parser | PASS；两个 Python 与统一入口无语法错误 | 可用 | 当前终端记录 |
| 05 正式 LOCO | `04_code/powershell/run_m06_pathway_loco.ps1` | EXIT 0；三次 rotation、Hallmark/Reactome 均完成 | 可用 | 本 run `results/`、`logs/` |
| 06 独立审计 | `M06_02_audit_pathway_loco.py` | PASS 24/24；no-leakage、输入、计数、core、边界声明全部通过 | 可用 | `results/M06_AUDIT_REPORT.md` |
| 07 图形 QC | 人工检查 6 张正式 PNG | PASS 6/6；无裁切、错位或不可读标签 | 可用 | `VISUAL_QC_REVIEW.PASS.txt` |
| 08 方法/结果文本冻结 | `M06_METHODS_RESULTS_LOCKED_v1.md` | PASS；方法、计数、限制、允许/禁止表述已写明 | 可用 | 本 run 根目录 |

## 下一步

复制到 `06_locked_results/modules/M06_LEAVE_ONE_COHORT_OUT/v1/`，更新模块登记与项目总进度；随后进入 M09 reporting freeze。M07/M08 为非 MVM 的 value-add，不阻塞 M09。
