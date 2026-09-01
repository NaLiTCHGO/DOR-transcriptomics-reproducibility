from __future__ import annotations

import argparse
import csv
import os
import sys
import zipfile
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def add_tree(zf: zipfile.ZipFile, root: Path, source: Path, prefix: str, manifest: list[dict[str, object]], excludes: tuple[str, ...] = ()) -> None:
    for file in sorted(source.rglob("*")):
        if not file.is_file():
            continue
        rel = file.relative_to(source).as_posix()
        if any(rel == x or rel.startswith(x.rstrip("/") + "/") for x in excludes):
            continue
        arc = f"{prefix}/{rel}" if rel else prefix
        zf.write(file, arc)
        category = prefix.split("/", 1)[-1].split("/", 1)[0]
        manifest.append({"package_path": arc, "source_path": file.relative_to(root).as_posix(), "bytes": file.stat().st_size, "category": category})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--workbook", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--zip", required=True)
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    workbook = Path(args.workbook).resolve()
    out = Path(args.output_dir).resolve()
    zip_path = Path(args.zip).resolve()
    out.mkdir(parents=True, exist_ok=True)
    zip_path.parent.mkdir(parents=True, exist_ok=True)

    modules = read_csv(root / "02_protocol_and_design/MODULE_REGISTRY.csv")
    scores = read_csv(root / "01_feasibility_and_decisions/decision_update_20260813/UPDATED_SCORECARD.csv")
    cohorts = read_csv(root / "06_locked_results/modules/M01_PROVENANCE_LOCK/v1_CORE_GEO/results/CORE_RNASEQ_COHORT_SUMMARY.csv")
    core = read_csv(root / "06_locked_results/modules/M06_LEAVE_ONE_COHORT_OUT/v1/results/M06_UNIVERSAL_LOCO_CORE_PATHWAYS.csv")
    loco = read_csv(root / "06_locked_results/modules/M06_LEAVE_ONE_COHORT_OUT/v1/results/M06_COLLECTION_LOCO_SUMMARY.csv")

    claims = [
        {"claim_id":"CLM_DECISION","class":"PERMITTED","claim":"The project remains GO_FULL under the frozen cohort-first design.","source":"01_feasibility_and_decisions/decision_update_20260813/FINAL_GO_DECISION.md","figure_table":"Executive summary"},
        {"claim_id":"CLM_REPRO","class":"PERMITTED","claim":"Gene-level effects show weak cross-cohort reproducibility and high heterogeneity.","source":"06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/M04_METHODS_RESULTS_LOCKED_v1.md","figure_table":"Figure 3"},
        {"claim_id":"CLM_PATHWAY","class":"PERMITTED","claim":"Pathway convergence is cohort-dependent and not universal.","source":"06_locked_results/modules/M05_PATHWAY_CONVERGENCE/v1/M05_METHODS_RESULTS_LOCKED_v1.md","figure_table":"Figure 4"},
        {"claim_id":"CLM_LOCO","class":"PERMITTED","claim":"A narrow eight-pathway universal LOCO core exists across the same three cohorts.","source":"06_locked_results/modules/M06_LEAVE_ONE_COHORT_OUT/v1/M06_METHODS_RESULTS_LOCKED_v1.md","figure_table":"Figure 5"},
        {"claim_id":"CLM_NO_EXTERNAL_VALIDATION","class":"PROHIBITED","claim":"Externally or clinically validated biomarker/pathway.","source":"06_locked_results/modules/M06_LEAVE_ONE_COHORT_OUT/v1/M06_METHODS_RESULTS_LOCKED_v1.md","figure_table":"All"},
        {"claim_id":"CLM_NO_ACTIVATION","class":"PROHIBITED","claim":"Negative NES proves pathway inhibition or causal mechanism.","source":"06_locked_results/modules/M05_PATHWAY_CONVERGENCE/v1/M05_METHODS_RESULTS_LOCKED_v1.md","figure_table":"All"},
        {"claim_id":"CLM_NO_BIOMARKER","class":"PROHIBITED","claim":"FMNL1, PGAP1, meta-FDR genes, or LOCO pathways are validated biomarkers.","source":"06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/M04_METHODS_RESULTS_LOCKED_v1.md","figure_table":"All"},
    ]
    claim_path = out / "M09_CLAIM_TRACEABILITY.csv"
    write_csv(claim_path, claims, ["claim_id","class","claim","source","figure_table"])

    pv = next(float(x["weighted_score"]) for x in scores if x["dimension"] == "project_value_total")
    er = next(float(x["weighted_score"]) for x in scores if x["dimension"] == "evidence_readiness_total")
    hall_core = sum(x["collection"] == "HALLMARK" for x in core)
    react_core = sum(x["collection"] == "REACTOME" for x in core)
    md = f"""# DOR 纯生信项目最终总览

`FINAL_STATE: GO_FULL / M00-M06 COMPLETE / M09 REPORTING FREEZE`

生成日期：2026-08-15（Asia/Shanghai）

## 一、最终结论

- 唯一正式决策：`GO_FULL`。
- Project Value Score：{pv:.1f}/100。
- Evidence Readiness Score：{er:.1f}/100（项目展示四舍五入为 86.5）。
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
3. M06：三轮均只用两个保留队列选择、第三队列留出验证；无阈值重调。Universal core 为 {hall_core} Hallmark + {react_core} Reactome，主紧凑结果为负向富集的 `HALLMARK_P53_PATHWAY`。
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
"""
    md_path = out / "DOR_Final_Project_Overview_20260815.md"
    md_path.write_text(md, encoding="utf-8")

    checks: list[dict[str, object]] = []
    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "pass": passed, "detail": detail})

    analysis_modules = {m["module_id"]: m["state"] for m in modules}
    for mid in [f"M{i:02d}_" for i in range(7)]:
        matches = [(k, v) for k, v in analysis_modules.items() if k.startswith(mid)]
        check(f"{mid}complete", len(matches) == 1 and matches[0][1].startswith("COMPLETE"), str(matches))
    m09_state = analysis_modules.get("M09_REPORTING_FREEZE", "MISSING")
    check("M09_reporting_freeze_complete", m09_state.startswith("COMPLETE"), m09_state)
    check("decision_GO_FULL", "GO_FULL" in (root / "01_feasibility_and_decisions/decision_update_20260813/FINAL_GO_DECISION.md").read_text(encoding="utf-8"), "authoritative decision file")
    check("project_value_85_1", abs(pv - 85.1) < 1e-9, str(pv))
    check("evidence_readiness_86_45", abs(er - 86.45) < 1e-9, str(er))
    check("three_core_cohorts", len(cohorts) == 3 and sum(int(x["N_independent"]) for x in cohorts) == 30, f"rows={len(cohorts)} independent_N={sum(int(x['N_independent']) for x in cohorts)}")
    check("m06_core_1_7", len(core) == 8 and hall_core == 1 and react_core == 7, f"total={len(core)} Hallmark={hall_core} Reactome={react_core}")
    check("m06_rotations_6", len(loco) == 6, f"rows={len(loco)}")
    check("workbook_exists", workbook.exists() and workbook.stat().st_size > 10_000, f"bytes={workbook.stat().st_size if workbook.exists() else 0}")
    check("markdown_boundaries", "不能泛化为全局通路一致" in md and "外部/临床/前瞻性验证" in md, "required claim boundaries")

    manifest: list[dict[str, object]] = []
    package_root = "DOR_Final_Scientific_Package_20260815"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        direct = [
            (workbook, f"{package_root}/overview/{workbook.name}", "overview"),
            (md_path, f"{package_root}/overview/{md_path.name}", "overview"),
            (claim_path, f"{package_root}/overview/{claim_path.name}", "overview"),
            (root / "PROJECT_PROGRESS.md", f"{package_root}/PROJECT_PROGRESS.md", "control"),
            (root / "ANALYSIS_WORKFLOW.md", f"{package_root}/ANALYSIS_WORKFLOW.md", "control"),
            (root / "00_project_control/LOW_TOKEN_WORKING_RULES.md", f"{package_root}/00_project_control/LOW_TOKEN_WORKING_RULES.md", "00_project_control"),
        ]
        for src, arc, category in direct:
            zf.write(src, arc)
            manifest.append({"package_path":arc,"source_path":src.relative_to(root).as_posix() if src.is_relative_to(root) else str(src),"bytes":src.stat().st_size,"category":category})
        add_tree(zf, root, root / "00_project_control/state_handoff", f"{package_root}/00_project_control/state_handoff", manifest)
        add_tree(zf, root, root / "01_feasibility_and_decisions", f"{package_root}/01_feasibility_and_decisions", manifest)
        add_tree(zf, root, root / "02_protocol_and_design", f"{package_root}/02_protocol_and_design", manifest)
        add_tree(zf, root, root / "04_code", f"{package_root}/04_code", manifest, excludes=("javascript/node_modules", "python/__pycache__"))
        add_tree(zf, root, root / "06_locked_results", f"{package_root}/06_locked_results", manifest, excludes=("modules/M09_REPORTING_FREEZE",))
        add_tree(zf, root, root / "07_manuscript", f"{package_root}/07_manuscript", manifest)

    manifest_path = out / "M09_RELEASE_MANIFEST.csv"
    write_csv(manifest_path, manifest, ["package_path","source_path","bytes","category"])
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        forbidden = [x for x in names if "/03_data/" in x or "/05_analysis_steps/" in x or x.lower().endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))]
        check("zip_created", zip_path.exists() and zip_path.stat().st_size > 1_000_000, f"bytes={zip_path.stat().st_size}")
        check("zip_no_raw_or_transient", not forbidden, f"forbidden={forbidden[:5]}")
        check("zip_has_workbook", any(x.endswith(workbook.name) for x in names), f"entries={len(names)}")
        check("zip_has_code", any("/04_code/python/M06_01_run_pathway_loco.py" in x for x in names), "M06 formal code")
        check("zip_has_locked_m06", any("/06_locked_results/modules/M06_LEAVE_ONE_COHORT_OUT/v1/M06_METHODS_RESULTS_LOCKED_v1.md" in x for x in names), "M06 locked narrative")

    audit_csv = out / "M09_AUDIT_CHECKS.csv"
    write_csv(audit_csv, checks, ["check","pass","detail"])
    failed = [x for x in checks if not x["pass"]]
    report = "# M09 Final Release Audit\n\n"
    report += f"`STATUS: {'FAIL' if failed else 'PASS'}`\n\n"
    report += f"- Checks passed: {len(checks)-len(failed)}/{len(checks)}.\n- Checks failed: {len(failed)}.\n"
    report += f"- ZIP entries: {len(manifest)}.\n- ZIP size: {zip_path.stat().st_size/1024/1024:.2f} MB.\n- Routine hashes: not used under the local working policy.\n\n"
    report += "## Failed checks\n\n" + ("None.\n" if not failed else "\n".join(f"- {x['check']}: {x['detail']}" for x in failed) + "\n")
    (out / "M09_AUDIT_REPORT.md").write_text(report, encoding="utf-8")
    if failed:
        print(report)
        return 1
    print(f"M09 release audit PASS: {len(checks)}/{len(checks)} checks; {len(manifest)} ZIP entries; {zip_path.stat().st_size/1024/1024:.2f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
