param(
  [string]$ProjectRoot=(Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path,
  [string]$PythonExe='python',
  [string]$RscriptExe='Rscript'
)

$ErrorActionPreference='Stop'
$env:DOR_PROJECT_ROOT=$ProjectRoot

$ErrorActionPreference='Continue'
& $RscriptExe --vanilla (Join-Path $ProjectRoot '04_code\R\M03_01_run_gse274832_deseq2.R')
$RExit=$LASTEXITCODE
$ErrorActionPreference='Stop'
if($RExit -ne 0){throw "M03 DESeq2 analysis failed with exit code $RExit"}
& $PythonExe (Join-Path $ProjectRoot '04_code\python\M03_01_audit_gse274832_candidate_stability.py')
if($LASTEXITCODE -ne 0){throw 'M03 candidate audit failed'}

Write-Host 'PASS: GSE274832 M03 statistical analysis reproduced'
