import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const args = Object.fromEntries(process.argv.slice(2).map((x, i, a) => x.startsWith("--") ? [x.slice(2), a[i + 1]] : null).filter(Boolean));
const projectRoot = path.resolve(args.projectRoot || process.cwd());
const outputPath = path.resolve(args.output || path.join(projectRoot, "outputs", "DOR_M09_20260815", "DOR_Project_Final_Overview_20260815.xlsx"));
const previewDir = path.resolve(args.previewDir || path.join(path.dirname(outputPath), "previews"));

function parseCsv(text) {
  const rows = [];
  let row = [], field = "", quoted = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"' && text[i + 1] === '"') { field += '"'; i++; }
      else if (ch === '"') quoted = false;
      else field += ch;
    } else if (ch === '"') quoted = true;
    else if (ch === ',') { row.push(field); field = ""; }
    else if (ch === '\n') { row.push(field.replace(/\r$/, "")); rows.push(row); row = []; field = ""; }
    else field += ch;
  }
  if (field.length || row.length) { row.push(field.replace(/\r$/, "")); rows.push(row); }
  return rows.filter(r => r.some(v => v !== ""));
}

async function readCsv(rel) {
  const rows = parseCsv(await fs.readFile(path.join(projectRoot, rel), "utf8"));
  const headers = rows.shift();
  headers[0] = headers[0].replace(/^\uFEFF/, "");
  return rows.map(r => Object.fromEntries(headers.map((h, i) => [h, r[i] ?? ""])));
}

function n(v) { const x = Number(v); return Number.isFinite(x) ? x : v; }
function colLetter(count) { let s = ""; while (count > 0) { count--; s = String.fromCharCode(65 + count % 26) + s; count = Math.floor(count / 26); } return s; }

const colors = {navy:"#17365D", teal:"#0F6B78", cyan:"#DDEBF7", pale:"#E2F0D9", amber:"#FFF2CC", red:"#FCE4D6", gray:"#F2F2F2", line:"#D9E2F3", white:"#FFFFFF", ink:"#1F2937"};

function titleBand(sheet, title, subtitle, lastCol) {
  sheet.showGridLines = false;
  sheet.getRange(`A1:${lastCol}1`).merge();
  sheet.getRange("A1").values = [[title]];
  sheet.getRange(`A1:${lastCol}1`).format = {fill:colors.navy,font:{bold:true,color:colors.white,size:16},verticalAlignment:"center"};
  sheet.getRange(`A1:${lastCol}1`).format.rowHeight = 30;
  sheet.getRange(`A2:${lastCol}2`).merge();
  sheet.getRange("A2").values = [[subtitle]];
  sheet.getRange(`A2:${lastCol}2`).format = {fill:colors.cyan,font:{color:colors.ink,italic:true},wrapText:true,verticalAlignment:"center"};
  sheet.getRange(`A2:${lastCol}2`).format.rowHeight = 32;
}

function writeTable(sheet, title, subtitle, headers, rows, widths, tableName) {
  const lastCol = colLetter(headers.length);
  titleBand(sheet, title, subtitle, lastCol);
  sheet.getRange(`A4:${lastCol}4`).values = [headers];
  sheet.getRange(`A4:${lastCol}4`).format = {fill:colors.teal,font:{bold:true,color:colors.white},wrapText:true,verticalAlignment:"center",borders:{preset:"outside",style:"thin",color:colors.navy}};
  sheet.getRange(`A4:${lastCol}4`).format.rowHeight = 30;
  if (rows.length) {
    sheet.getRange(`A5:${lastCol}${rows.length + 4}`).values = rows;
    sheet.getRange(`A5:${lastCol}${rows.length + 4}`).format = {font:{color:colors.ink,size:10},wrapText:true,verticalAlignment:"top",borders:{insideHorizontal:{style:"thin",color:colors.line}}};
    sheet.getRange(`A5:${lastCol}${rows.length + 4}`).format.rowHeight = 42;
    const table = sheet.tables.add(`A4:${lastCol}${rows.length + 4}`, true, tableName);
    table.style = "TableStyleMedium2";
    table.showBandedRows = true;
  }
  widths.forEach((w, i) => { sheet.getRange(`${colLetter(i + 1)}:${colLetter(i + 1)}`).format.columnWidth = w; });
  sheet.freezePanes.freezeRows(4);
}

const modules = await readCsv("02_protocol_and_design/MODULE_REGISTRY.csv");
const scoresRaw = await readCsv("01_feasibility_and_decisions/decision_update_20260813/UPDATED_SCORECARD.csv");
const coreCohorts = await readCsv("06_locked_results/modules/M01_PROVENANCE_LOCK/v1_CORE_GEO/results/CORE_RNASEQ_COHORT_SUMMARY.csv");
const pairwise = await readCsv("06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/results/M04_PAIRWISE_REPRODUCIBILITY.csv");
const m04Counts = await readCsv("06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/results/M04_KEY_COUNTS.csv");
const pathwayPairs = await readCsv("06_locked_results/modules/M05_PATHWAY_CONVERGENCE/v1/results/M05_PATHWAY_PAIRWISE_REPRODUCIBILITY.csv");
const strictPaths = await readCsv("06_locked_results/modules/M05_PATHWAY_CONVERGENCE/v1/results/M05_STRICT_CONVERGENT_PATHWAYS.csv");
const loco = await readCsv("06_locked_results/modules/M06_LEAVE_ONE_COHORT_OUT/v1/results/M06_COLLECTION_LOCO_SUMMARY.csv");
const core = await readCsv("06_locked_results/modules/M06_LEAVE_ONE_COHORT_OUT/v1/results/M06_UNIVERSAL_LOCO_CORE_PATHWAYS.csv");
const reqs = await readCsv("02_protocol_and_design/REQUIREMENTS_TRACEABILITY_MATRIX.csv");
const debts = await readCsv("02_protocol_and_design/DEBT_REGISTER.csv");

const workbook = Workbook.create();
const names = ["Executive Summary","Module Status","Cohort Freeze","Scores","Cross-cohort","Pathways and LOCO","Claims and Limits","Traceability","Debt Register","README"];
const sheets = Object.fromEntries(names.map(name => [name, workbook.worksheets.add(name)]));

const moduleRows = modules.map(m => {
  const analysisComplete = /^M0[0-6]_/.test(m.module_id) && m.state.startsWith("COMPLETE") ? 1 : 0;
  return [m.module_id,m.name,m.criticality,m.MVM_required === "True",m.state,m.depends_on,m.contract_path,m.dod_path,analysisComplete];
});
writeTable(sheets["Module Status"],"Module status and critical path","M00–M06 are locked analysis modules. M07/M08 are optional value-add; M09 is the reporting-freeze module.",["Module","Name","Criticality","MVM required","State","Depends on","Contract","DoD / accepted marker","M00–M06 complete flag"],moduleRows,[20,34,14,14,30,24,42,44,20],"ModuleStatusTable");

const cohortRows = coreCohorts.map(c => [c.series,c.N_recruited_context,n(c.N_omics),n(c.N_independent),n(c.DOR_omics),n(c.NOR_omics),"Human ovarian granulosa cells","Core RNA-seq",c.m01_verdict,c.series === "GSE274832" ? "Retina/mm10/mm8 text conflict resolved as metadata-template error; disclose." : "Cohort-specific QC/model limitations remain."]);
cohortRows.push(["E-MTAB-391","13 DOR + 13 NOR participants reported; 28 cycle-samples",28,26,14,14,"Granulosa cells","Legacy sensitivity","PASS_WITH_LIMITATION","28 cycle-samples from 26 patients; repeated-patient Source IDs unavailable."]);
writeTable(sheets["Cohort Freeze"],"Frozen accession and three-N table","N_recruited/context, N_omics, and N_independent are kept separate. DOR, POI, and POR labels are not interchanged.",["Accession","N recruited / context","N omics","N independent","DOR omics","NOR omics","Tissue","Role","Verdict","Provenance / limitation"],cohortRows,[18,34,12,14,12,12,24,20,22,54],"CohortFreezeTable");

const valueTotal = scoresRaw.find(x => x.dimension === "project_value_total");
const readinessTotal = scoresRaw.find(x => x.dimension === "evidence_readiness_total");
const scoreRows = [
  ["value","Project Value Score",n(valueTotal.weight),n(valueTotal.score_0_to_10),n(valueTotal.weighted_score),valueTotal.update_basis],
  ["readiness","Evidence Readiness Score",n(readinessTotal.weight),n(readinessTotal.score_0_to_10),n(readinessTotal.weighted_score),readinessTotal.update_basis],
  ...scoresRaw.filter(x => !["project_value_total","evidence_readiness_total"].includes(x.dimension)).map(x => [x.score_family,x.dimension,n(x.weight),n(x.score_0_to_10),n(x.weighted_score),x.update_basis])
];
writeTable(sheets["Scores"],"Final feasibility scorecard","The two top rows are the authoritative project-facing totals; component rows retain the final delta reassessment basis.",["Family","Dimension","Weight","Score (0–10)","Weighted score (/100)","Update basis"],scoreRows,[14,42,12,14,20,68],"ScorecardTable");
sheets["Scores"].getRange(`C5:E${scoreRows.length+4}`).format.numberFormat = "0.00";

const crossRows = [
  ...pairwise.map(x => ["M04","Pairwise gene effect",`${x.cohort_a} vs ${x.cohort_b}`,n(x.spearman_rho_log2fc),n(x.direction_concordance_all_nonzero),"Weak gene-level reproducibility","06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/results/M04_PAIRWISE_REPRODUCIBILITY.csv"]),
  ["M04","Common gene universe","All three cohorts",13993,"","Frozen intersection","06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/results/M04_KEY_COUNTS.csv"],
  ["M04","Random-effects FDR<0.05","All three cohorts",449,"","Exploratory; not biomarkers","06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/results/M04_KEY_COUNTS.csv"],
  ["M04","I2>=75%","All three cohorts",5883,"","High heterogeneity; k=3 descriptive","06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/results/M04_KEY_COUNTS.csv"],
  ["M04","Strict exploratory consensus","FMNL1 / PGAP1",2,"","Exploratory effects only","06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/M04_METHODS_RESULTS_LOCKED_v1.md"]
];
writeTable(sheets["Cross-cohort"],"Cross-cohort reproducibility and heterogeneity","All effects are DOR−NOR and were estimated within cohort before synthesis; no megamatrix pooling or clinical classifier was used.",["Module","Metric","Comparison / object","Value","Rate","Interpretation","Locked source"],crossRows,[12,28,32,14,14,34,72],"CrossCohortTable");
sheets["Cross-cohort"].getRange(`D5:E${crossRows.length+4}`).format.numberFormat = "0.0000";
sheets["Cross-cohort"].getRange("E5:E7").format.numberFormat = "0.0%";
sheets["Cross-cohort"].getRange("D8:D11").format.numberFormat = "#,##0";

const pathwayRows = [
  ...pathwayPairs.map(x => ["M05",x.collection,`${x.cohort_a} vs ${x.cohort_b}`,"Pathway NES rho",n(x.pathway_nes_spearman_rho),n(x.pathway_direction_concordance),"Pairwise pathway context"]),
  ["M05","HALLMARK","All three","Strict convergent pathways",8,"","All negative-in-DOR rank"],
  ["M05","REACTOME","All three","Strict convergent pathways",47,"","44 negative, 3 positive; terms are redundant"],
  ...loco.map(x => ["M06",x.collection,`hold out ${x.held_out_cohort}`,"selected → direction → strict",`${x.retained_pair_selected_n} → ${x.held_out_direction_replication_n} → ${x.held_out_strict_replication_n}`,n(x.strict_replication_rate_among_selected),"Held-out cohort was not used for selection"]),
  ...core.map(x => ["M06",x.collection,"Universal LOCO core",x.pathway,n(x.median_abs_nes),3,"Negative-in-DOR rank; internal conditional robustness"])
];
writeTable(sheets["Pathways and LOCO"],"Pathway convergence and leave-one-cohort-out robustness","MSigDB Human v2026.1.Hs Hallmark is primary; Reactome is secondary. Negative/positive NES is enrichment direction, not pathway inhibition/activation.",["Module","Collection","Comparison","Metric / pathway","Value","Rate / rotations","Interpretation boundary"],pathwayRows,[12,14,26,52,18,18,54],"PathwayLocoTable");
sheets["Pathways and LOCO"].getRange("E5:E10").format.numberFormat = "0.0000";
sheets["Pathways and LOCO"].getRange("E11:E12").format.numberFormat = "#,##0";
sheets["Pathways and LOCO"].getRange("F13:F18").format.numberFormat = "0.0%";
sheets["Pathways and LOCO"].getRange("E19:E26").format.numberFormat = "0.000";
sheets["Pathways and LOCO"].getRange("F19:F26").format.numberFormat = "0";

const claims = [
  ["PERMITTED","Project decision","The project remains GO_FULL within the frozen cohort-first reproducibility design.","Project Value 85.1; Evidence Readiness 86.5."],
  ["PERMITTED","Gene-level result","Gene-level effects show weak cross-cohort reproducibility and high heterogeneity.","Report pairwise rho and I2 with k=3 limitations."],
  ["PERMITTED","Pathway result","A narrow universal LOCO pathway core exists across the same three cohorts.","Primary compact core: HALLMARK_P53_PATHWAY."],
  ["PERMITTED","Cohort limitation","GSE232306 is the most discordant held-out cohort for Hallmark replication.","Describe phenotype-aligned technical/latent confounding."],
  ["PROHIBITED","Validation","Externally/clinically/prospectively validated biomarker or pathway.","M06 is internal conditional robustness only."],
  ["PROHIBITED","Mechanism","Negative NES proves pathway inhibition or causality.","NES is rank enrichment direction."],
  ["PROHIBITED","Biomarker","FMNL1, PGAP1, 449 meta-FDR genes, or 8 core pathways are validated biomarkers.","They remain exploratory signals."],
  ["PROHIBITED","Generality","The three cohorts show global pathway concordance.","Most M05 strict pathways survive only one retained-pair rotation."],
  ["PROHIBITED","Phenotype","DOR, POI, and POR are interchangeable.","Preserve original phenotype definitions and provenance."]
];
writeTable(sheets["Claims and Limits"],"Frozen claims and interpretation boundaries","These statements govern project summaries, figures, manuscripts, and reuse of the locked results.",["Class","Topic","Frozen statement","Required qualifier"],claims,[16,24,72,62],"ClaimsTable");
sheets["Claims and Limits"].getRange(`A5:A${claims.length+4}`).conditionalFormats.add("containsText",{text:"PERMITTED",format:{fill:colors.pale,font:{bold:true,color:"#375623"}}});
sheets["Claims and Limits"].getRange(`A5:A${claims.length+4}`).conditionalFormats.add("containsText",{text:"PROHIBITED",format:{fill:colors.red,font:{bold:true,color:"#9C0006"}}});

const traceRows = reqs.map(r => [r.requirement_id,r.requirement_type,r.requirement_text,r.module_id,r.result_artifact_id,r.figure_table_id,r.manuscript_claim_id,r.lock_id,r.status]);
writeTable(sheets["Traceability"],"Requirement-to-result traceability","Every scientific and release requirement points to a module, locked artifact, intended figure/table, claim, and lock level.",["Requirement","Type","Requirement text","Module","Result artifact","Figure / table","Claim","Lock","Status"],traceRows,[20,14,62,24,28,22,24,18,28],"TraceabilityTable");

const debtRows = debts.map(d => [d.debt_id,d.type,d.severity,d.origin,d.affected_object,d.description,d.owner,d.due_before_gate,d.blocking === "True",d.status,d.closure_evidence]);
writeTable(sheets["Debt Register"],"Scientific and engineering debt register","OPEN items are non-blocking but must remain visible. CONTROLLED means bounded by design or reporting, not eliminated.",["Debt","Type","Severity","Origin","Affected object","Description","Owner","Due gate","Blocking","Status","Closure / control evidence"],debtRows,[12,14,12,28,30,62,16,14,12,16,72],"DebtTable");
sheets["Debt Register"].getRange(`J5:J${debtRows.length+4}`).conditionalFormats.add("containsText",{text:"OPEN",format:{fill:colors.amber,font:{bold:true,color:"#7F6000"}}});
sheets["Debt Register"].getRange(`J5:J${debtRows.length+4}`).conditionalFormats.add("containsText",{text:"CLOSED",format:{fill:colors.pale,font:{bold:true,color:"#375623"}}});

const readmeRows = [
  ["Workbook role","Final project overview and traceability index; it does not replace locked CSVs or methods files."],
  ["Scientific question","Phenotype/provenance-aware cross-cohort reproducibility and heterogeneity of DOR transcriptomic programs."],
  ["Build date","2026-08-15"],
  ["Decision authority","01_feasibility_and_decisions/decision_update_20260813/FINAL_GO_DECISION.md"],
  ["Daily recovery entry","PROJECT_PROGRESS.md"],
  ["Locked result root","06_locked_results/modules/"],
  ["Formal code root","04_code/"],
  ["No routine hashes","Per local single-user working policy; versioned folders, logs, audits, and locked paths are used."],
  ["Primary gene-set source","https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/"],
  ["MSigDB collection page","https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp"],
  ["External-validation boundary","No large same-tissue independent clinical validation cohort is included."],
  ["Reproducibility boundary","The delivery ZIP excludes raw FASTQ and transient run directories; download/access manifests remain in the project."],
];
writeTable(sheets["README"],"Workbook README","Read this sheet first. Numeric source-of-truth remains the cited locked artifact.",["Item","Value / path"],readmeRows,[30,110],"ReadmeTable");

const dash = sheets["Executive Summary"];
dash.showGridLines = false;
dash.getRange("A1:H1").merge(); dash.getRange("A1").values = [["DOR transcriptomics reproducibility project — final overview"]];
dash.getRange("A1:H1").format = {fill:colors.navy,font:{bold:true,color:colors.white,size:18},verticalAlignment:"center"}; dash.getRange("A1:H1").format.rowHeight=34;
dash.getRange("A2:H2").merge(); dash.getRange("A2").values = [["Final reporting freeze after M00–M06; generated 2026-08-15. Read with the locked-result limitations."]];
dash.getRange("A2:H2").format = {fill:colors.cyan,font:{italic:true,color:colors.ink},verticalAlignment:"center"}; dash.getRange("A2:H2").format.rowHeight=26;
dash.getRange("A4:B4").merge(); dash.getRange("A4").values=[["Final decision"]];
dash.getRange("C4:D4").merge(); dash.getRange("C4").values=[["Project Value"]];
dash.getRange("E4:F4").merge(); dash.getRange("E4").values=[["Evidence Readiness"]];
dash.getRange("G4:H4").merge(); dash.getRange("G4").values=[["Completed analysis modules"]];
dash.getRange("A5:B6").merge(); dash.getRange("A5").values=[["GO_FULL"]];
dash.getRange("C5:D6").merge(); dash.getRange("C5").formulas=[["='Scores'!E5"]];
dash.getRange("E5:F6").merge(); dash.getRange("E5").formulas=[["='Scores'!E6"]];
dash.getRange("G5:H6").merge(); dash.getRange("G5").formulas=[["=SUM('Module Status'!$I$5:$I$15)"]];
dash.getRange("A4:H4").format={fill:colors.teal,font:{bold:true,color:colors.white},horizontalAlignment:"center"};
dash.getRange("A5:H6").format={fill:colors.pale,font:{bold:true,color:colors.navy,size:18},horizontalAlignment:"center",verticalAlignment:"center",borders:{preset:"outside",style:"thin",color:colors.teal}};
dash.getRange("C5:F6").format.numberFormat="0.0";
dash.getRange("A8:H8").merge(); dash.getRange("A8").values=[["What the completed analysis supports"]]; dash.getRange("A8:H8").format={fill:colors.navy,font:{bold:true,color:colors.white}};
const highlights=[
  ["Data foundation","Three core granulosa-cell RNA-seq cohorts: 6 + 12 + 12 independent omics samples; E-MTAB-391 remains legacy sensitivity only."],
  ["Gene level","13,993 common genes; pairwise effect rho 0.0491–0.2363; 5,883 genes have descriptive I2>=75%."],
  ["Pathway level","GSE274832/GSE193136 show better pathway than gene-rank concordance, but this is not universal across GSE232306 pairs."],
  ["Robust core","Universal LOCO core contains 1 Hallmark + 7 Reactome pathways; compact primary result is negative enrichment of P53 PATHWAY."],
  ["Decision","GO_FULL remains valid because provenance, QC, cohort effects, synthesis, pathway analysis, and internal robustness gates all closed."],
];
dash.getRange("A9:B13").merge(true); dash.getRange("A9:A13").values=highlights.map(x=>[x[0]]);
dash.getRange("C9:H13").merge(true); dash.getRange("C9:C13").values=highlights.map(x=>[x[1]]);
dash.getRange("A9:B13").format={fill:colors.gray,font:{bold:true,color:colors.navy},verticalAlignment:"top",wrapText:true};
dash.getRange("C9:H13").format={font:{color:colors.ink},verticalAlignment:"top",wrapText:true,borders:{insideHorizontal:{style:"thin",color:colors.line}}};
dash.getRange("A15:H15").merge(); dash.getRange("A15").values=[["Non-negotiable interpretation boundary"]]; dash.getRange("A15:H15").format={fill:"#9C0006",font:{bold:true,color:colors.white}};
dash.getRange("A16:H18").merge(); dash.getRange("A16").values=[["This project does not provide external clinical validation, diagnostic performance, causal pathway proof, or evidence that negative NES equals functional inhibition. FMNL1, PGAP1, meta-FDR genes, and the eight-pathway LOCO core remain exploratory scientific signals."]];
dash.getRange("A16:H18").format={fill:colors.red,font:{color:"#9C0006",bold:true},wrapText:true,verticalAlignment:"center"};
dash.getRange("A20:H20").merge(); dash.getRange("A20").values=[["Next publication step: use the Traceability and Claims and Limits sheets to assemble figures/tables and refresh the literature collision audit immediately before submission."]];
dash.getRange("A20:H20").format={fill:colors.amber,font:{color:"#7F6000",italic:true},wrapText:true}; dash.getRange("A20:H20").format.rowHeight=36;
for (let i=1;i<=8;i++) dash.getRange(`${colLetter(i)}:${colLetter(i)}`).format.columnWidth=18;
dash.freezePanes.freezeRows(2);

await fs.mkdir(path.dirname(outputPath), {recursive:true});
await fs.mkdir(previewDir, {recursive:true});
const keyInspect = await workbook.inspect({kind:"table",range:"Executive Summary!A1:H20",include:"values,formulas",tableMaxRows:20,tableMaxCols:8,maxChars:5000});
console.log(keyInspect.ndjson);
const errors = await workbook.inspect({kind:"match",searchTerm:"#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",options:{useRegex:true,maxResults:300},summary:"final formula error scan",maxChars:3000});
console.log(errors.ndjson);
for (const name of names) {
  const blob = await workbook.render({sheetName:name,autoCrop:"all",scale:1,format:"png"});
  await fs.writeFile(path.join(previewDir, `${name.replaceAll(" ","_")}.png`), new Uint8Array(await blob.arrayBuffer()));
}
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
console.log(JSON.stringify({status:"PASS",outputPath,previewDir,sheets:names.length,strictPathways:strictPaths.length,m04KeyCountRows:m04Counts.length},null,2));
