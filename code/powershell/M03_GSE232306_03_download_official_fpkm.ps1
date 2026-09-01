param([string]$ProjectRoot=$env:DOR_PROJECT_ROOT)
$ErrorActionPreference='Stop'
if(-not $ProjectRoot){$ProjectRoot=(Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path}

$OutDir=Join-Path $ProjectRoot '03_data\processed_external_READ_ONLY\GSE232306'
$OutFile=Join-Path $OutDir 'GSE232306_1_genes_fpkm_expression.xlsx'
$Url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232306/suppl/GSE232306_1_genes_fpkm_expression.xlsx'
New-Item -ItemType Directory -Force -Path $OutDir|Out-Null
if(-not(Test-Path -LiteralPath $OutFile)){
  & curl.exe -L --fail --retry 3 --output $OutFile $Url
  if($LASTEXITCODE-ne0){throw "Official FPKM download failed with exit $LASTEXITCODE"}
}
if((Get-Item -LiteralPath $OutFile).Length-lt 10000){throw 'Official FPKM workbook is unexpectedly small'}
Write-Output "PASS: official GSE232306 FPKM workbook available at $OutFile"
