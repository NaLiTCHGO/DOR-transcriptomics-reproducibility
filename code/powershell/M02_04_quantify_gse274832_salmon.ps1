param([string]$ProjectRoot=$env:DOR_PROJECT_ROOT)
$ErrorActionPreference='Stop'
if(-not $ProjectRoot){$ProjectRoot=(Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path}
$RunRoot=Join-Path $ProjectRoot '05_analysis_steps\M02_COHORT_QC\runs\20260813_M02_B1_GSE274832'
$Manifest=@(Import-Csv (Join-Path $ProjectRoot '04_code\configs\GSE274832_FASTQ_MANIFEST.csv'))
$FastqRoot=Join-Path $ProjectRoot '03_data\raw_external_READ_ONLY\GSE274832\fastq'
foreach($File in $Manifest){$File|Add-Member -NotePropertyName local_path -NotePropertyValue (Join-Path $FastqRoot $File.filename) -Force}
$Locked=@(Import-Csv (Join-Path $ProjectRoot '04_code\configs\GSE274832_SAMPLES.csv'))
$Index=Join-Path $ProjectRoot '03_data\references\GENCODE_v50_GRCh38p14\salmon_index_v1.11.4'
$OutRoot=Join-Path $RunRoot 'results\salmon'
$LogRoot=Join-Path $RunRoot 'logs\salmon'
$Salmon=Join-Path $ProjectRoot '00_project_control\tools\linux_x86_64\salmon-1.11.4\bin\salmon'
New-Item -ItemType Directory -Force -Path $OutRoot,$LogRoot|Out-Null
if(-not(Test-Path (Join-Path $Index 'versionInfo.json'))){throw 'Validated Salmon index is not ready'}

function To-Wsl([string]$Path){'/mnt/'+$Path.Substring(0,1).ToLower()+'/'+$Path.Substring(3).Replace('\','/')}
$Summary=@()
foreach($Sample in $Locked){
    $Srr=$Sample.srr
    $R1=($Manifest|Where-Object {$_.srr-eq$Srr -and $_.mate-eq'1'}|Select-Object -First 1).local_path
    $R2=($Manifest|Where-Object {$_.srr-eq$Srr -and $_.mate-eq'2'}|Select-Object -First 1).local_path
    if(-not(Test-Path $R1) -or -not(Test-Path $R2)){throw "Missing FASTQ pair: $Srr"}
    $Out=Join-Path $OutRoot $Srr
    if(-not(Test-Path (Join-Path $Out 'quant.sf'))){
        if(Test-Path $Out){
            $FailedOut="$Out.failed_"+(Get-Date -Format yyyyMMdd_HHmmss)
            Move-Item $Out $FailedOut
        }
        New-Item -ItemType Directory -Force -Path $Out|Out-Null
        $Args=@('--exec',(To-Wsl $Salmon),'quant','-i',(To-Wsl $Index),'-l','A','-1',(To-Wsl $R1),'-2',(To-Wsl $R2),'-p','8','--validateMappings','--gcBias','--seqBias','-o',(To-Wsl $Out))
        # Salmon checks its online version endpoint at startup.  An offline
        # "Version Server Response: Not Found" message is written to stderr
        # but is not a quantification failure, so judge the native exit code.
        $ErrorActionPreference='Continue'
        & wsl.exe @Args 1>> (Join-Path $LogRoot "$Srr.stdout.log") 2>> (Join-Path $LogRoot "$Srr.stderr.log")
        $SalmonExit=$LASTEXITCODE
        $ErrorActionPreference='Stop'
        if($SalmonExit-ne 0){throw "Salmon quant failed: $Srr (exit $SalmonExit)"}
    }
    $MetaPath=Join-Path $Out 'aux_info\meta_info.json'
    if(-not(Test-Path $MetaPath)){throw "Missing Salmon meta_info: $Srr"}
    $Meta=Get-Content $MetaPath -Raw|ConvertFrom-Json
    $Pct=[double]$Meta.percent_mapped
    $Verdict=if($Pct-ge 60){'PASS'}elseif($Pct-ge 40){'REVIEW'}else{'FAIL'}
    $Summary+=[pscustomobject]@{
        srr=$Srr;phenotype=$Sample.phenotype;gsm=$Sample.gsm
        num_processed=[long]$Meta.num_processed;num_mapped=[long]$Meta.num_mapped
        percent_mapped=$Pct;library_types=($Meta.library_types-join ';')
        mapping_qc_verdict=$Verdict
    }
    $Summary|Export-Csv (Join-Path $OutRoot 'SALMON_MAPPING_QC.csv') -NoTypeInformation -Encoding UTF8
}
if(($Summary|Where-Object {$_.mapping_qc_verdict-eq'FAIL'}).Count-gt 0){throw 'One or more samples failed the mapping-rate screen'}
$Summary|Sort-Object phenotype,srr|Format-Table -AutoSize
Write-Output 'PASS: 6/6 GSE274832 samples quantified; mapping summary written'
