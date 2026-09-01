param(
  [string]$ProjectRoot=(Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path,
  [string]$PythonExe='python',
  [string]$RscriptExe='C:\Program Files\R\R-4.5.3\bin\Rscript.exe'
)
$ErrorActionPreference='Stop'
$env:DOR_PROJECT_ROOT=$ProjectRoot
$M02=Join-Path $ProjectRoot '05_analysis_steps\M02_COHORT_QC\runs\20260814_M02_B3_GSE193136'
$M03=Join-Path $ProjectRoot '05_analysis_steps\M03_WITHIN_COHORT_EFFECTS\runs\20260814_M03_B3_GSE193136'
New-Item -ItemType Directory -Force -Path (Join-Path $M02 'logs'),(Join-Path $M02 'results'),(Join-Path $M03 'logs'),(Join-Path $M03 'results')|Out-Null
$Transcript=Join-Path $M02 'logs\GSE193136_PIPELINE_TRANSCRIPT.log'
Start-Transcript -Path $Transcript -Append|Out-Null

function Stage([string]$Name,[scriptblock]$Action){
  Write-Host "START $Name $(Get-Date -Format o)"
  & $Action
  if($LASTEXITCODE -ne 0){throw "$Name failed with exit code $LASTEXITCODE"}
  Write-Host "PASS  $Name $(Get-Date -Format o)"
}

try {
  if(Test-Path (Join-Path $M02 'STEP01_FASTQ_INTEGRITY.PASS.txt')){Write-Host 'SKIP  01 FASTQ integrity: PASS marker present'}else{Stage '01 FASTQ expected-size and gzip integrity' {& (Join-Path $ProjectRoot '04_code\powershell\M02_GSE193136_01_validate_fastq.ps1') -ProjectRoot $ProjectRoot}}
  if(Test-Path (Join-Path $M02 'STEP02_FASTP.PASS.txt')){Write-Host 'SKIP  02 fastp: PASS marker present'}else{Stage '02 fastp raw read QC' {& (Join-Path $ProjectRoot '04_code\powershell\M02_GSE193136_02_run_fastp.ps1') -ProjectRoot $ProjectRoot}}
  if(Test-Path (Join-Path $M02 'STEP03_RAW_QC.PASS.txt')){Write-Host 'SKIP  03 raw-QC summary: PASS marker present'}else{Stage '03 summarize raw QC' {& $PythonExe (Join-Path $ProjectRoot '04_code\python\M02_GSE193136_03_summarize_raw_qc.py')}}
  if(Test-Path (Join-Path $M02 'STEP04_SALMON.PASS.txt')){Write-Host 'SKIP  04 Salmon: PASS marker present'}else{Stage '04 Salmon quantification' {& (Join-Path $ProjectRoot '04_code\powershell\M02_GSE193136_04_run_salmon.ps1') -ProjectRoot $ProjectRoot}}
  if(Test-Path (Join-Path $M02 'STEP05_EXPRESSION_QC.PASS.txt')){Write-Host 'SKIP  05 expression QC: PASS marker present'}else{Stage '05 expression matrices and numerical QC' {& $PythonExe (Join-Path $ProjectRoot '04_code\python\M02_GSE193136_05_build_expression_qc.py')}}
  if(Test-Path (Join-Path $M02 'STEP06_QC_PLOTS.PASS.txt')){Write-Host 'SKIP  06 QC plots: PASS marker present'}else{Stage '06 expression QC plots' {& $RscriptExe --vanilla (Join-Path $ProjectRoot '04_code\R\M02_GSE193136_06_plot_expression_qc.R')}}
  if(Test-Path (Join-Path $M02 'STEP07_M02_FINAL_GATE.PASS.txt')){Write-Host 'SKIP  07 M02 final gate: PASS marker present'}else{Stage '07 M02 final gate' {& $PythonExe (Join-Path $ProjectRoot '04_code\python\M02_GSE193136_07_finalize_qc.py')}}
  if(Test-Path (Join-Path $M03 'GSE193136_M03_PASS_WITH_LIMITATION.txt')){Write-Host 'SKIP  08 M03 DESeq2: PASS marker present'}else{Stage '08 M03 age-adjusted DESeq2 and LOSO' {& $RscriptExe --vanilla (Join-Path $ProjectRoot '04_code\R\M03_GSE193136_01_run_age_adjusted_deseq2.R')}}
  if(Test-Path (Join-Path $M03 'STEP09_M03_CANDIDATE_AUDIT.PASS.txt')){Write-Host 'SKIP  09 candidate audit: PASS marker present'}else{Stage '09 M03 candidate stability audit' {& $PythonExe (Join-Path $ProjectRoot '04_code\python\M03_GSE193136_02_audit_candidate_stability.py')}}
  Write-Host 'PASS: GSE193136 M02-to-M03 workflow completed'
} finally {
  Stop-Transcript|Out-Null
}
