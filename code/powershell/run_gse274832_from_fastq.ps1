param(
  [string]$ProjectRoot=(Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path,
  [string]$PythonExe='python',
  [string]$RscriptExe='Rscript',
  [string]$ProxyUrl=$env:DOR_HTTP_PROXY
)

$ErrorActionPreference='Stop'
$env:DOR_PROJECT_ROOT=$ProjectRoot

function Stage([string]$Name,[scriptblock]$Action){
  Write-Host "START $Name"
  & $Action
  if($LASTEXITCODE -ne 0){throw "$Name failed with exit code $LASTEXITCODE"}
  Write-Host "PASS  $Name"
}

Stage 'FASTQ integrity and fastp' {
  & (Join-Path $ProjectRoot '04_code\powershell\M02_01_validate_fastq_and_run_fastp.ps1') -ProjectRoot $ProjectRoot
}
Stage 'Raw-QC summary' {
  & (Join-Path $ProjectRoot '04_code\powershell\M02_02_summarize_gse274832_raw_qc.ps1') -ProjectRoot $ProjectRoot
}
Stage 'GENCODE reference and Salmon index' {
  & (Join-Path $ProjectRoot '04_code\powershell\M02_03_prepare_gencode50_reference.ps1') -ProjectRoot $ProjectRoot -ProxyUrl $ProxyUrl
}
Stage 'Salmon quantification' {
  & (Join-Path $ProjectRoot '04_code\powershell\M02_04_quantify_gse274832_salmon.ps1') -ProjectRoot $ProjectRoot
}
Stage 'Gene matrices and numerical QC' {
  & $PythonExe (Join-Path $ProjectRoot '04_code\python\M02_01_build_gse274832_expression_qc.py')
}
Stage 'QC plots' {
  $ErrorActionPreference='Continue'
  & $RscriptExe --vanilla (Join-Path $ProjectRoot '04_code\R\M02_01_plot_gse274832_expression_qc.R')
  $script:RNativeExit=$LASTEXITCODE
  $ErrorActionPreference='Stop'
  if($script:RNativeExit -ne 0){throw "QC plots failed with exit code $script:RNativeExit"}
}
Stage 'M02 final gate' {
  & $PythonExe (Join-Path $ProjectRoot '04_code\python\M02_02_finalize_gse274832_qc.py')
}
Stage 'M03 DESeq2 and leave-one-sample-out analysis' {
  $ErrorActionPreference='Continue'
  & $RscriptExe --vanilla (Join-Path $ProjectRoot '04_code\R\M03_01_run_gse274832_deseq2.R')
  $script:RNativeExit=$LASTEXITCODE
  $ErrorActionPreference='Stop'
  if($script:RNativeExit -ne 0){throw "M03 DESeq2 failed with exit code $script:RNativeExit"}
}
Stage 'M03 FDR-candidate audit' {
  & $PythonExe (Join-Path $ProjectRoot '04_code\python\M03_01_audit_gse274832_candidate_stability.py')
}

Write-Host 'PASS: GSE274832 reproducibility workflow completed from FASTQ'
