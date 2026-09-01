import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const root = process.env.DOR_PROJECT_ROOT;
if (!root) throw new Error("DOR_PROJECT_ROOT is required");
const input = path.join(root, "03_data/processed_external_READ_ONLY/GSE232306/GSE232306_1_genes_fpkm_expression.xlsx");
const outDir = path.join(root, "05_analysis_steps/M03_WITHIN_COHORT_EFFECTS/runs/20260814_M03_B4_GSE232306/results/official_fpkm_validation");
await fs.mkdir(outDir, { recursive: true });
const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(input));
const summary = await workbook.inspect({kind: "workbook,sheet,table", maxChars: 8000, tableMaxRows: 20, tableMaxCols: 20, tableMaxCellChars: 100});
const sheets = await workbook.inspect({kind: "sheet", include: "id,name", maxChars: 3000});
await fs.writeFile(path.join(outDir, "GSE232306_OFFICIAL_FPKM_WORKBOOK_INSPECT.ndjson"), `${summary.ndjson}\n${sheets.ndjson}\n`, "utf8");
process.stdout.write(`${summary.ndjson}\n${sheets.ndjson}\n`);
