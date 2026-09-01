param([string]$ProjectRoot=$env:DOR_PROJECT_ROOT)
$ErrorActionPreference='Stop'
if(-not $ProjectRoot){$ProjectRoot=(Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path}
function To-Wsl([string]$Path){'/mnt/'+$Path.Substring(0,1).ToLower()+'/'+$Path.Substring(3).Replace('\','/')}

$RunRoot=Join-Path $ProjectRoot '05_analysis_steps\M02_COHORT_QC\runs\20260814_M02_B4_GSE232306'
$Manifest=@(Import-Csv (Join-Path $ProjectRoot '04_code\configs\GSE232306_FASTQ_MANIFEST.csv'))
$Samples=@(Import-Csv (Join-Path $ProjectRoot '04_code\configs\GSE232306_SAMPLES.csv'))
$FastqRoot=Join-Path $ProjectRoot '03_data\raw_external_READ_ONLY\GSE232306\fastq'
$Index=Join-Path $ProjectRoot '03_data\references\GENCODE_v50_GRCh38p14\salmon_index_v1.11.4'
$Salmon=Join-Path $ProjectRoot '00_project_control\tools\linux_x86_64\salmon-1.11.4\bin\salmon'
$OutRoot=Join-Path $RunRoot 'results\salmon'
$LogRoot=Join-Path $RunRoot 'logs\salmon'
New-Item -ItemType Directory -Force -Path $OutRoot,$LogRoot | Out-Null
if(-not(Test-Path (Join-Path $RunRoot 'STEP03_RAW_QC.PASS.txt'))){throw 'STEP03 raw-QC marker missing'}
if(-not(Test-Path (Join-Path $Index 'versionInfo.json'))){throw 'Validated Salmon index is missing'}

$Summary=@()
foreach($Sample in $Samples){
    $Srr=$Sample.srr
    $R1=Join-Path $FastqRoot (($Manifest|Where-Object {$_.srr-eq$Srr -and $_.mate-eq'1'}|Select-Object -First 1).filename)
    $R2=Join-Path $FastqRoot (($Manifest|Where-Object {$_.srr-eq$Srr -and $_.mate-eq'2'}|Select-Object -First 1).filename)
    $Out=Join-Path $OutRoot $Srr
    if(-not(Test-Path (Join-Path $Out 'quant.sf'))){
        if(Test-Path $Out){Move-Item -LiteralPath $Out -Destination ($Out+'.failed_'+(Get-Date -Format yyyyMMdd_HHmmss))}
        New-Item -ItemType Directory -Force -Path $Out|Out-Null
        $Args=@('--exec',(To-Wsl $Salmon),'quant','-i',(To-Wsl $Index),'-l','A','-1',(To-Wsl $R1),'-2',(To-Wsl $R2),'-p','8','--validateMappings','--gcBias','--seqBias','-o',(To-Wsl $Out))
        $ErrorActionPreference='Continue'
        & wsl.exe @Args 1>> (Join-Path $LogRoot "$Srr.stdout.log") 2>> (Join-Path $LogRoot "$Srr.stderr.log")
        $NativeExit=$LASTEXITCODE
        $ErrorActionPreference='Stop'
        if($NativeExit-ne 0){throw "Salmon failed: $Srr (exit $NativeExit)"}
    }
    $MetaPath=Join-Path $Out 'aux_info\meta_info.json'
    if(-not(Test-Path $MetaPath)){throw "Missing Salmon meta_info: $Srr"}
    $Meta=Get-Content $MetaPath -Raw|ConvertFrom-Json
    $Pct=[double]$Meta.percent_mapped
    # Decoy-aware transcriptome mapping >=60% passes the mapping screen.
    # Rates from 30%-60% remain REVIEW and require cohort-wide expression
    # coherence and granulosa-identity QC; <30% is a hard screen failure.
    $Verdict=if($Pct-ge60){'PASS'}elseif($Pct-ge30){'REVIEW_LOW_TRANSCRIPTOME_MAPPING'}else{'FAIL'}
    $Summary += [pscustomobject]@{srr=$Srr;gsm=$Sample.gsm;phenotype=$Sample.phenotype;num_processed=[long]$Meta.num_processed;num_mapped=[long]$Meta.num_mapped;percent_mapped=$Pct;library_types=($Meta.library_types-join';');mapping_qc_verdict=$Verdict}
    $Summary|Export-Csv (Join-Path $OutRoot 'GSE232306_SALMON_MAPPING_QC.csv') -NoTypeInformation -Encoding UTF8
}
if($Summary.Count-ne12){throw "Expected 12 quantifications; found $($Summary.Count)"}
if(($Summary|Where-Object mapping_qc_verdict -eq 'FAIL').Count-gt0){throw 'One or more samples had <30% transcriptome mapping and failed the hard screen'}
"PASS`t12/12 Salmon quantifications`t$(Get-Date -Format o)" | Set-Content (Join-Path $RunRoot 'STEP04_SALMON.PASS.txt') -Encoding UTF8
$Summary|Sort-Object phenotype,srr|Format-Table -AutoSize
Write-Output 'PASS: GSE232306 Salmon quantification completed for 12/12 samples'
