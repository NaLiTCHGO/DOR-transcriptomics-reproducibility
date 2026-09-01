param(
  [string]$ProjectRoot=(Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path,
  [string]$PythonExe='python',
  [string]$RscriptExe='Rscript'
)

$ErrorActionPreference='Stop'
$env:DOR_PROJECT_ROOT=$ProjectRoot

$required=Get-ChildItem (Join-Path $ProjectRoot '05_analysis_steps\M02_COHORT_QC\runs\20260813_M02_B1_GSE274832\results\salmon') -Recurse -Filter quant.sf -File -ErrorAction SilentlyContinue
if(@($required).Count -ne 6){throw "Expected 6 Salmon quant.sf files; found $(@($required).Count)"}

& $PythonExe (Join-Path $ProjectRoot '04_code\python\M02_01_build_gse274832_expression_qc.py')
if($LASTEXITCODE -ne 0){throw 'Expression-matrix reconstruction failed'}
$ErrorActionPreference='Continue'
& $RscriptExe --vanilla (Join-Path $ProjectRoot '04_code\R\M02_01_plot_gse274832_expression_qc.R')
$RExit=$LASTEXITCODE
$ErrorActionPreference='Stop'
if($RExit -ne 0){throw "QC plotting failed with exit code $RExit"}
& $PythonExe (Join-Path $ProjectRoot '04_code\python\M02_02_finalize_gse274832_qc.py')
if($LASTEXITCODE -ne 0){throw 'M02 final gate failed'}
$ErrorActionPreference='Continue'
& $RscriptExe --vanilla (Join-Path $ProjectRoot '04_code\R\M03_01_run_gse274832_deseq2.R')
$RExit=$LASTEXITCODE
$ErrorActionPreference='Stop'
if($RExit -ne 0){throw "M03 DESeq2 analysis failed with exit code $RExit"}
& $PythonExe (Join-Path $ProjectRoot '04_code\python\M03_01_audit_gse274832_candidate_stability.py')
if($LASTEXITCODE -ne 0){throw 'M03 candidate audit failed'}

Write-Host 'PASS: GSE274832 workflow reproduced from existing Salmon estimates'
