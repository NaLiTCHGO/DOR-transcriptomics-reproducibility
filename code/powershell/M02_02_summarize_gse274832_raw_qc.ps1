param([string]$ProjectRoot=$env:DOR_PROJECT_ROOT)
$ErrorActionPreference='Stop'
if(-not $ProjectRoot){$ProjectRoot=(Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path}
$RunRoot=Join-Path $ProjectRoot '05_analysis_steps\M02_COHORT_QC\runs\20260813_M02_B1_GSE274832'
$Input=Import-Csv (Join-Path $RunRoot 'results\fastp\FASTP_RAW_QC_ALL6_SUMMARY.csv')
$Locked=Import-Csv (Join-Path $ProjectRoot '04_code\configs\GSE274832_SAMPLES.csv')
$Out=Join-Path $RunRoot 'results\cohort_raw_qc'; New-Item -ItemType Directory -Force -Path $Out|Out-Null
foreach($x in $Input){$CurrentSrr=$x.srr;$x.phenotype=($Locked|Where-Object {$_.srr -eq $CurrentSrr}|Select-Object -First 1).phenotype}
$Input|Export-Csv (Join-Path $RunRoot 'results\fastp\FASTP_RAW_QC_ALL6_SUMMARY.csv') -NoTypeInformation -Encoding UTF8
$Input|Export-Csv (Join-Path $Out 'GSE274832_FASTP_SAMPLE_QC.csv') -NoTypeInformation -Encoding UTF8
$Group=$Input|Group-Object phenotype|ForEach-Object{[pscustomobject]@{phenotype=$_.Name;n=$_.Count;mean_reads=[math]::Round((($_.Group|Measure-Object total_reads -Average).Average),0);mean_q30=[math]::Round((($_.Group|Measure-Object q30_rate -Average).Average),4);min_q30=[math]::Round((($_.Group|Measure-Object q30_rate -Minimum).Minimum),4);mean_gc=[math]::Round((($_.Group|Measure-Object gc_content -Average).Average),4);mean_duplication=[math]::Round((($_.Group|Measure-Object duplication_rate -Average).Average),4)}}
$Group|Export-Csv (Join-Path $Out 'GSE274832_FASTP_GROUP_SUMMARY.csv') -NoTypeInformation -Encoding UTF8
$Qc=$Input|ForEach-Object{[pscustomobject]@{srr=$_.srr;phenotype=$_.phenotype;q30_rate=[double]$_.q30_rate;gc_content=[double]$_.gc_content;duplication_rate=[double]$_.duplication_rate;raw_qc_verdict=if([double]$_.q30_rate-ge 0.85 -and [double]$_.gc_content-ge 0.35 -and [double]$_.gc_content-le 0.65){'PASS'}else{'REVIEW'}}}
$Qc|Export-Csv (Join-Path $Out 'GSE274832_RAW_QC_VERDICT.csv') -NoTypeInformation -Encoding UTF8
if(($Qc|Where-Object raw_qc_verdict-ne'PASS').Count-gt 0){throw 'One or more samples require raw QC review'}
Write-Output 'PASS: 6/6 GSE274832 samples passed prespecified raw-QC screen'
$Group|Format-Table -AutoSize
