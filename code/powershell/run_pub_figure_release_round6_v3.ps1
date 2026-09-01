param(
  [string]$ProjectRoot = ".",
  [string]$SourceReleaseRunId = "20260824_RBE_FIG2_S6_V03",
  [string]$ReleaseRunId = "20260824_RBE_FIG2_S6_V04",
  [string]$SkillRoot = "external\research_scientific_figure_suite",
  [string]$Rscript = "C:\Program Files\R\R-4.5.3\bin\Rscript.exe"
)

$ErrorActionPreference = "Stop"
$ProjectRoot = (Resolve-Path -LiteralPath $ProjectRoot).Path
$SkillRoot = (Resolve-Path -LiteralPath $SkillRoot).Path
$sourceRun = Join-Path $ProjectRoot "05_analysis_steps\PUB_FIGURE_RELEASE\runs\$SourceReleaseRunId"
$releaseRun = Join-Path $ProjectRoot "05_analysis_steps\PUB_FIGURE_RELEASE\runs\$ReleaseRunId"
foreach ($required in @($sourceRun, $Rscript)) {
  if (-not (Test-Path -LiteralPath $required)) { throw "Missing required input: $required" }
}
if (Test-Path -LiteralPath $releaseRun) { throw "Release run already exists: $releaseRun" }

foreach ($folder in @(
    $releaseRun,
    (Join-Path $releaseRun "source_data"),
    (Join-Path $releaseRun "reconstructed_heatmaps"),
    (Join-Path $releaseRun "figures\main"),
    (Join-Path $releaseRun "figures\supplement"),
    (Join-Path $releaseRun "qc"),
    (Join-Path $releaseRun "manifests"),
    (Join-Path $releaseRun "logs"),
    (Join-Path $releaseRun "code_snapshot")
  )) {
  New-Item -ItemType Directory -Path $folder -Force | Out-Null
}

Get-ChildItem -LiteralPath (Join-Path $sourceRun "source_data") | ForEach-Object {
  Copy-Item -LiteralPath $_.FullName -Destination (Join-Path $releaseRun "source_data") -Recurse
}
Get-ChildItem -LiteralPath (Join-Path $sourceRun "reconstructed_heatmaps") | ForEach-Object {
  Copy-Item -LiteralPath $_.FullName -Destination (Join-Path $releaseRun "reconstructed_heatmaps") -Recurse
}
foreach ($relative in @(
    "qc\VST_HEATMAP_RECONSTRUCTION_VALIDATION.csv",
    "qc\RUNTIME_WARNING_LOG_VST_RECONSTRUCTION.csv",
    "qc\RUNTIME_WARNING_LOG_SUPP_COHORTS.csv",
    "qc\FULL_RUNTIME_WARNING_RECONCILIATION_SUPP_COHORTS.md",
    "manifests\VST_HEATMAP_RECONSTRUCTION_PROVENANCE.md",
    "manifests\R_SESSION_INFO_VST_RECONSTRUCTION.txt"
  )) {
  $src = Join-Path $sourceRun $relative
  if (-not (Test-Path -LiteralPath $src)) { throw "Missing carried QC/provenance file: $src" }
  Copy-Item -LiteralPath $src -Destination (Join-Path $releaseRun (Split-Path $relative -Parent))
}

# Only Figure 2 and Figure S4 are revised in V3. The other eight figures are
# copied byte-identically from the locked V2 parent release.
$carriedFigureIds = @("Figure_3", "Figure_4", "Figure_5", "Figure_S1", "Figure_S2", "Figure_S3", "Figure_S5", "Figure_S6")
foreach ($figureId in $carriedFigureIds) {
  $category = if ($figureId -match '^Figure_[3-5]$') { "main" } else { "supplement" }
  Copy-Item -LiteralPath (Join-Path $sourceRun "figures\$category\$figureId") -Destination (Join-Path $releaseRun "figures\$category") -Recurse
}

$carryRows = foreach ($figureId in $carriedFigureIds) {
  $category = if ($figureId -match '^Figure_[3-5]$') { "main" } else { "supplement" }
  $sourceFolder = Join-Path $sourceRun "figures\$category\$figureId"
  $releaseFolder = Join-Path $releaseRun "figures\$category\$figureId"
  Get-ChildItem -LiteralPath $releaseFolder -File | ForEach-Object {
    $sourceFile = Join-Path $sourceFolder $_.Name
    $sourceHash = (Get-FileHash -LiteralPath $sourceFile -Algorithm SHA256).Hash.ToLowerInvariant()
    $releaseHash = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    [pscustomobject]@{
      figure_id = $figureId
      artifact = $_.Name
      source_release_run = $SourceReleaseRunId
      bytes = $_.Length
      source_sha256 = $sourceHash
      release_sha256 = $releaseHash
      status = if ($sourceHash -eq $releaseHash) { "CARRIED_FORWARD_BYTE_IDENTICAL" } else { "FAIL_HASH_MISMATCH" }
    }
  }
}
$carryRows | Export-Csv -LiteralPath (Join-Path $releaseRun "manifests\CARRIED_FORWARD_ARTIFACTS.csv") -NoTypeInformation -Encoding UTF8
if ($carryRows | Where-Object status -ne "CARRIED_FORWARD_BYTE_IDENTICAL") {
  throw "One or more carried figure artifacts differ from the V2 parent release"
}

$scripts = @(
  "04_code\R\PUB_FIG2_02_render_candidate.R",
  "04_code\R\PUB_FIGBATCH_00_common.R",
  "04_code\R\PUB_FIGBATCH_04_render_supp_synthesis.R",
  "04_code\python\PUB_RELEASE_01_validate_cross_format.py",
  "04_code\python\PUB_RELEASE_02_validate_legends.py",
  "04_code\powershell\run_pub_figure_release_qc_v1.ps1",
  "04_code\powershell\run_pub_figure_release_round6_v3.ps1"
)
foreach ($relative in $scripts) {
  Copy-Item -LiteralPath (Join-Path $ProjectRoot $relative) -Destination (Join-Path $releaseRun "code_snapshot")
}

function Invoke-RStep {
  param([string]$Label, [string]$Script, [string[]]$Arguments)
  $log = Join-Path $releaseRun ("logs\" + $Label + ".log")
  & $Rscript --vanilla $Script @Arguments 2>&1 | Tee-Object -FilePath $log
  if ($LASTEXITCODE -ne 0) { throw "R step failed: $Label" }
}

Invoke-RStep -Label "01_render_figure2_release_v3" -Script (Join-Path $ProjectRoot "04_code\R\PUB_FIG2_02_render_candidate.R") -Arguments @(
  "--project-root", $ProjectRoot, "--run-dir", $releaseRun, "--skill-root", $SkillRoot, "--release-mode", "TRUE"
)
Invoke-RStep -Label "02_render_s4_release_v3" -Script (Join-Path $ProjectRoot "04_code\R\PUB_FIGBATCH_04_render_supp_synthesis.R") -Arguments @(
  "--project-root=$ProjectRoot", "--run-root=$releaseRun", "--skill-root=$SkillRoot", "--figure-set=S4", "--release-mode=TRUE"
)

@(
  "# RBE result-figure release revision 3",
  "",
  "- Run ID: ``$ReleaseRunId``",
  "- Parent release: ``$SourceReleaseRunId``",
  "- Re-rendered: Figure 2 (shared A-F header baseline) and Figure S4 (right-side label 9 and common-universe subtitle).",
  "- Carried forward byte-identically: Figure 3-Figure 5, Figure S1-S3, and Figure S5-S6.",
  "- Figure S4 count verification: 5,717 = 5,698 + 10 + 9 within the 13,993-gene common universe.",
  "- Scientific-result recomputation: NO",
  "- User candidate visual decision: PASS",
  "- Status: RENDER_COMPLETE_QC_PENDING"
) | Set-Content -LiteralPath (Join-Path $releaseRun "RUN_STATUS.md") -Encoding UTF8

Write-Output "PASS: revision-3 final figure bundle rendered at $releaseRun"
