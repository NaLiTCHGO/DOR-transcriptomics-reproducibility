param(
  [string]$ProjectRoot = ".",
  [string]$SourceReleaseRunId = "20260824_RBE_FIG2_S6_V03",
  [string]$CandidateRunId = "20260824_VISUAL_QC_ROUND6_HEADER_S4",
  [string]$SkillRoot = "external\research_scientific_figure_suite",
  [string]$Rscript = "C:\Program Files\R\R-4.5.3\bin\Rscript.exe"
)

$ErrorActionPreference = "Stop"
$ProjectRoot = (Resolve-Path -LiteralPath $ProjectRoot).Path
$SkillRoot = (Resolve-Path -LiteralPath $SkillRoot).Path
$sourceRun = Join-Path $ProjectRoot "05_analysis_steps\PUB_FIGURE_RELEASE\runs\$SourceReleaseRunId"
$candidateRun = Join-Path $ProjectRoot "05_analysis_steps\PUB_FIGURE_REBUILD\runs\$CandidateRunId"
foreach ($required in @($sourceRun, $Rscript)) {
  if (-not (Test-Path -LiteralPath $required)) { throw "Missing required input: $required" }
}
if (Test-Path -LiteralPath $candidateRun) { throw "Candidate run already exists: $candidateRun" }

foreach ($folder in @(
    $candidateRun,
    (Join-Path $candidateRun "source_data"),
    (Join-Path $candidateRun "reconstructed_heatmaps"),
    (Join-Path $candidateRun "figures"),
    (Join-Path $candidateRun "qc"),
    (Join-Path $candidateRun "manifests"),
    (Join-Path $candidateRun "logs"),
    (Join-Path $candidateRun "code_snapshot")
  )) {
  New-Item -ItemType Directory -Path $folder -Force | Out-Null
}
Get-ChildItem -LiteralPath (Join-Path $sourceRun "source_data") | ForEach-Object {
  Copy-Item -LiteralPath $_.FullName -Destination (Join-Path $candidateRun "source_data") -Recurse
}
Get-ChildItem -LiteralPath (Join-Path $sourceRun "reconstructed_heatmaps") | ForEach-Object {
  Copy-Item -LiteralPath $_.FullName -Destination (Join-Path $candidateRun "reconstructed_heatmaps") -Recurse
}

$scripts = @(
  "04_code\R\PUB_FIG2_02_render_candidate.R",
  "04_code\R\PUB_FIGBATCH_00_common.R",
  "04_code\R\PUB_FIGBATCH_04_render_supp_synthesis.R",
  "04_code\powershell\run_pub_visual_qc_round6_microrepair.ps1"
)
foreach ($relative in $scripts) {
  Copy-Item -LiteralPath (Join-Path $ProjectRoot $relative) -Destination (Join-Path $candidateRun "code_snapshot")
}

function Invoke-RStep {
  param([string]$Label, [string]$Script, [string[]]$Arguments)
  $log = Join-Path $candidateRun ("logs\" + $Label + ".log")
  & $Rscript --vanilla $Script @Arguments 2>&1 | Tee-Object -FilePath $log
  if ($LASTEXITCODE -ne 0) { throw "R step failed: $Label" }
}

Invoke-RStep -Label "01_render_figure2_header_alignment" -Script (Join-Path $ProjectRoot "04_code\R\PUB_FIG2_02_render_candidate.R") -Arguments @(
  "--project-root", $ProjectRoot, "--run-dir", $candidateRun, "--skill-root", $SkillRoot, "--release-mode", "FALSE"
)
Invoke-RStep -Label "02_render_s4_label_and_scope_repair" -Script (Join-Path $ProjectRoot "04_code\R\PUB_FIGBATCH_04_render_supp_synthesis.R") -Arguments @(
  "--project-root=$ProjectRoot", "--run-root=$candidateRun", "--skill-root=$SkillRoot", "--figure-set=S4", "--release-mode=FALSE"
)

@(
  "# Visual QC round 6 microrepair",
  "",
  "- Source: locked result-figure release ``$SourceReleaseRunId``",
  "- Figure 2 change: rebuild every A-F header as one shared grob so the tag and accession use identical font metrics and a common vertical anchor.",
  "- Figure S4 change: move the count 9 horizontally to the right of its point without vertical displacement.",
  "- Figure S4 clarification: state that overlap counts use the 13,993-gene common universe.",
  "- Count verification: GSE232306 has 5,717 candidates in that universe = 5,698 unique + 10 shared with GSE274832 + 9 shared with GSE193136.",
  "- Scientific-result recomputation: NO",
  "- Candidate formats: 300-ppi QC PNG only",
  "- Status: CANDIDATE_READY_FOR_ACTUAL_RENDER_QC"
) | Set-Content -LiteralPath (Join-Path $candidateRun "RUN_ALIGNMENT.md") -Encoding UTF8

Write-Output "PASS: round-6 Figure 2/S4 candidates rendered at $candidateRun"
