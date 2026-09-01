param(
  [string]$ProjectRoot = ".",
  [string]$ReleaseRunId = "20260824_RBE_FIG2_S6_V01",
  [string]$Python = "python",
  [string]$PdfToPpm = "pdftoppm"
)

$ErrorActionPreference = "Stop"
$ProjectRoot = (Resolve-Path -LiteralPath $ProjectRoot).Path
$runRoot = Join-Path $ProjectRoot "05_analysis_steps\PUB_FIGURE_RELEASE\runs\$ReleaseRunId"
$renderRoot = Join-Path $runRoot "qc\pdf_rendered_300dpi"
foreach ($required in @($runRoot, $Python, $PdfToPpm)) {
  if (-not (Test-Path -LiteralPath $required)) { throw "Missing required input: $required" }
}
New-Item -ItemType Directory -Path $renderRoot -Force | Out-Null

$validator = Join-Path $ProjectRoot "04_code\python\PUB_RELEASE_01_validate_cross_format.py"
$thisScript = Join-Path $ProjectRoot "04_code\powershell\run_pub_figure_release_qc_v1.ps1"
foreach ($script in @($validator, $thisScript)) {
  Copy-Item -LiteralPath $script -Destination (Join-Path $runRoot "code_snapshot") -Force
}

$pdfs = Get-ChildItem -LiteralPath (Join-Path $runRoot "figures") -Recurse -Filter "*_submission.pdf" | Sort-Object FullName
if ($pdfs.Count -ne 10) { throw "Expected 10 final PDFs, found $($pdfs.Count)" }
foreach ($pdf in $pdfs) {
  $figureId = $pdf.Directory.Name
  $prefix = Join-Path $renderRoot ($figureId + "_pdf_300dpi")
  & $PdfToPpm -png -r 300 -singlefile $pdf.FullName $prefix
  if ($LASTEXITCODE -ne 0) { throw "PDF render failed: $($pdf.FullName)" }
}

& $Python $validator --project-root $ProjectRoot --run-root $runRoot 2>&1 |
  Tee-Object -FilePath (Join-Path $runRoot "logs\06_cross_format_qc.log")
if ($LASTEXITCODE -ne 0) { throw "Cross-format QC failed; see the QC table and log" }

Write-Output "PASS: all final PDFs opened and rendered; raster, geometry, font and format gates passed"
