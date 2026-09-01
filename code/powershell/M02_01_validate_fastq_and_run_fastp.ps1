param([string]$ProjectRoot=$env:DOR_PROJECT_ROOT)
$ErrorActionPreference='Stop'
if(-not $ProjectRoot){$ProjectRoot=(Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path}
function To-Wsl([string]$Path){'/mnt/'+$Path.Substring(0,1).ToLower()+'/'+$Path.Substring(3).Replace('\','/')}
$RunRoot=Join-Path $ProjectRoot '05_analysis_steps\M02_COHORT_QC\runs\20260813_M02_B1_GSE274832'
$Manifest=Import-Csv (Join-Path $ProjectRoot '04_code\configs\GSE274832_FASTQ_MANIFEST.csv')
$FastqRoot=Join-Path $ProjectRoot '03_data\raw_external_READ_ONLY\GSE274832\fastq'
foreach($File in $Manifest){$File|Add-Member -NotePropertyName local_path -NotePropertyValue (Join-Path $FastqRoot $File.filename) -Force}
$ResultRoot=Join-Path $RunRoot 'results\fastp'; $LogRoot=Join-Path $RunRoot 'logs\fastp'
$Fastp=To-Wsl (Join-Path $ProjectRoot '00_project_control\tools\linux_x86_64\fastp-1.3.6')
$Integrity=@()
foreach($Group in ($Manifest|Group-Object srr)){
    $Json=Join-Path $ResultRoot ($Group.Name+'.fastp.json')
    $Pair=$Group.Group|Sort-Object {[int]$_.mate}
    if($Pair.Count-ne 2){throw "Expected paired FASTQ files for $($Group.Name)"}
    foreach($File in $Pair){
        if(-not(Test-Path $File.local_path)){throw "Missing $($File.local_path)"}
        $Bytes=(Get-Item $File.local_path).Length; $Md5=(Get-FileHash $File.local_path -Algorithm MD5).Hash.ToLowerInvariant()
        if($Bytes-ne[long]$File.expected_bytes -or $Md5-ne$File.expected_md5.ToLowerInvariant()){throw "Integrity failed $($File.local_path)"}
        $Wsl=To-Wsl $File.local_path; wsl.exe --exec gzip -t $Wsl
        if($LASTEXITCODE-ne 0){throw "gzip failed $($File.local_path)"}
        $Integrity += [pscustomobject]@{srr=$Group.Name;phenotype=$File.phenotype;mate=$File.mate;file=$File.local_path;bytes=$Bytes;md5=$Md5;integrity='PASS'}
    }
    if(Test-Path $Json){continue}
    $R1=To-Wsl $Pair[0].local_path; $R2=To-Wsl $Pair[1].local_path
    $Html=Join-Path $ResultRoot ($Group.Name+'.fastp.html'); $Json=Join-Path $ResultRoot ($Group.Name+'.fastp.json'); $Log=Join-Path $LogRoot ($Group.Name+'.fastp.log')
    $WH=To-Wsl $Html; $WJ=To-Wsl $Json; $WL=To-Wsl $Log
    wsl.exe --exec bash -lc "'$Fastp' -i '$R1' -I '$R2' --json '$WJ' --html '$WH' --report_title '$($Group.Name) raw FASTQ QC' --thread 8 --disable_adapter_trimming > '$WL' 2>&1"
    if($LASTEXITCODE-ne 0){throw "fastp failed $($Group.Name)"}
}
$Integrity|Sort-Object srr,{[int]$_.mate}|Export-Csv (Join-Path $ResultRoot 'FASTQ_FILE_INTEGRITY_ALL6.csv') -NoTypeInformation -Encoding UTF8
$Summary=@(); foreach($Jf in Get-ChildItem $ResultRoot -Filter '*.fastp.json'){$J=Get-Content $Jf.FullName -Raw|ConvertFrom-Json;$Srr=$Jf.BaseName-replace '\.fastp$','';$Ph=($Manifest|Where-Object {$_.srr -eq $Srr}|Select-Object -First 1).phenotype;$Summary += [pscustomobject]@{srr=$Srr;phenotype=$Ph;total_reads=[long]$J.summary.before_filtering.total_reads;total_bases=[long]$J.summary.before_filtering.total_bases;q20_rate=[double]$J.summary.before_filtering.q20_rate;q30_rate=[double]$J.summary.before_filtering.q30_rate;gc_content=[double]$J.summary.before_filtering.gc_content;read1_mean_length=[double]$J.summary.before_filtering.read1_mean_length;read2_mean_length=[double]$J.summary.before_filtering.read2_mean_length;duplication_rate=[double]$J.duplication.rate;fastp_status='PASS'}}
$Summary|Sort-Object srr|Export-Csv (Join-Path $ResultRoot 'FASTP_RAW_QC_ALL6_SUMMARY.csv') -NoTypeInformation -Encoding UTF8
if($Summary.Count-ne 6 -or $Integrity.Count-ne 12){throw "Expected 6 reports/12 files; got $($Summary.Count)/$($Integrity.Count)"}
Write-Output 'PASS: 6/6 samples and 12/12 files complete'
$Summary|Format-Table -AutoSize
