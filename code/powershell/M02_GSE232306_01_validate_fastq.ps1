param([string]$ProjectRoot=$env:DOR_PROJECT_ROOT)
$ErrorActionPreference='Stop'
if(-not $ProjectRoot){$ProjectRoot=(Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path}

function To-Wsl([string]$Path){'/mnt/'+$Path.Substring(0,1).ToLower()+'/'+$Path.Substring(3).Replace('\','/')}

$RunRoot=Join-Path $ProjectRoot '05_analysis_steps\M02_COHORT_QC\runs\20260814_M02_B4_GSE232306'
$Manifest=@(Import-Csv (Join-Path $ProjectRoot '04_code\configs\GSE232306_FASTQ_MANIFEST.csv'))
$Samples=@(Import-Csv (Join-Path $ProjectRoot '04_code\configs\GSE232306_SAMPLES.csv'))
$FastqRoot=Join-Path $ProjectRoot '03_data\raw_external_READ_ONLY\GSE232306\fastq'
$OutRoot=Join-Path $RunRoot 'results\fastq_integrity'
$LogRoot=Join-Path $RunRoot 'logs\integrity'
New-Item -ItemType Directory -Force -Path $OutRoot,$LogRoot | Out-Null

if($Samples.Count-ne 12){throw "Expected 12 locked samples; found $($Samples.Count)"}
if(($Samples|Where-Object {$_.phenotype-eq'DOR'}).Count-ne 6 -or ($Samples|Where-Object {$_.phenotype-eq'NOR'}).Count-ne 6){throw 'Expected 6 DOR and 6 NOR samples'}
if($Manifest.Count-ne 24){throw "Expected 24 manifest rows; found $($Manifest.Count)"}

$Rows=@()
foreach($Group in ($Manifest|Group-Object srr)){
    $Pair=@($Group.Group|Sort-Object {[int]$_.mate})
    if($Pair.Count-ne 2 -or $Pair[0].mate-ne'1' -or $Pair[1].mate-ne'2'){throw "Invalid paired-end manifest for $($Group.Name)"}
    foreach($File in $Pair){
        $Path=Join-Path $FastqRoot $File.filename
        if(-not(Test-Path -LiteralPath $Path)){throw "Missing FASTQ: $Path"}
        $Bytes=(Get-Item -LiteralPath $Path).Length
        if($Bytes-ne[long]$File.expected_bytes){throw "File-size mismatch: $($File.filename); expected $($File.expected_bytes), observed $Bytes"}
        & wsl.exe --exec gzip -t (To-Wsl $Path) 2>> (Join-Path $LogRoot 'gzip_test.stderr.log')
        if($LASTEXITCODE-ne 0){throw "gzip integrity failed: $($File.filename)"}
        $Rows += [pscustomobject]@{srr=$File.srr;gsm=$File.gsm;phenotype=$File.phenotype;mate=[int]$File.mate;filename=$File.filename;expected_bytes=[long]$File.expected_bytes;observed_bytes=$Bytes;size_match='PASS';gzip_integrity='PASS';overall='PASS'}
    }
}
$Rows|Sort-Object srr,mate|Export-Csv (Join-Path $OutRoot 'GSE232306_FASTQ_INTEGRITY.csv') -NoTypeInformation -Encoding UTF8
"PASS`t$($Rows.Count)/24 FASTQ; size and gzip checks passed`t$(Get-Date -Format o)" | Set-Content (Join-Path $RunRoot 'STEP01_FASTQ_INTEGRITY.PASS.txt') -Encoding UTF8
Write-Output 'PASS: GSE232306 24/24 FASTQ passed expected-size and gzip integrity checks (MD5 intentionally omitted by local policy)'

