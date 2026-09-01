param([string]$ProjectRoot=$env:DOR_PROJECT_ROOT)
$ErrorActionPreference='Stop'
if(-not $ProjectRoot){$ProjectRoot=(Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path}
function To-Wsl([string]$Path){'/mnt/'+$Path.Substring(0,1).ToLower()+'/'+$Path.Substring(3).Replace('\','/')}

$RunRoot=Join-Path $ProjectRoot '05_analysis_steps\M02_COHORT_QC\runs\20260814_M02_B3_GSE193136'
$Manifest=@(Import-Csv (Join-Path $ProjectRoot '04_code\configs\GSE193136_FASTQ_MANIFEST.csv'))
$FastqRoot=Join-Path $ProjectRoot '03_data\raw_external_READ_ONLY\GSE193136\fastq'
$OutRoot=Join-Path $RunRoot 'results\fastp'
$LogRoot=Join-Path $RunRoot 'logs\fastp'
$Fastp=Join-Path $ProjectRoot '00_project_control\tools\linux_x86_64\fastp-1.3.6'
New-Item -ItemType Directory -Force -Path $OutRoot,$LogRoot | Out-Null
if(-not(Test-Path (Join-Path $RunRoot 'STEP01_FASTQ_INTEGRITY.PASS.txt'))){throw 'STEP01 integrity marker missing'}

foreach($Group in ($Manifest|Group-Object srr)){
    $Pair=@($Group.Group|Sort-Object {[int]$_.mate})
    $Json=Join-Path $OutRoot ($Group.Name+'.fastp.json')
    $Html=Join-Path $OutRoot ($Group.Name+'.fastp.html')
    $StdoutLog=Join-Path $LogRoot ($Group.Name+'.fastp.stdout.log')
    $StderrLog=Join-Path $LogRoot ($Group.Name+'.fastp.stderr.log')
    if(Test-Path $Json){
        try {
            $Existing=Get-Content -LiteralPath $Json -Raw|ConvertFrom-Json
            if($Existing.summary.before_filtering.total_reads){continue}
            throw 'Incomplete fastp JSON'
        } catch {
            Move-Item -LiteralPath $Json -Destination ($Json+'.failed_'+(Get-Date -Format yyyyMMdd_HHmmss))
            if(Test-Path $Html){Move-Item -LiteralPath $Html -Destination ($Html+'.failed_'+(Get-Date -Format yyyyMMdd_HHmmss))}
        }
    }
    $R1=Join-Path $FastqRoot $Pair[0].filename
    $R2=Join-Path $FastqRoot $Pair[1].filename
    $Args=@('--exec',(To-Wsl $Fastp),'-i',(To-Wsl $R1),'-I',(To-Wsl $R2),'--json',(To-Wsl $Json),'--html',(To-Wsl $Html),'--report_title',"$($Group.Name) raw FASTQ QC",'--thread','8','--disable_adapter_trimming')
    # fastp writes its normal progress and summary to stderr. PowerShell 5.1
    # can promote native stderr records to terminating errors when the global
    # preference is Stop, so capture the native exit code explicitly.
    $ErrorActionPreference='Continue'
    & wsl.exe @Args 1>> $StdoutLog 2>> $StderrLog
    $FastpExit=$LASTEXITCODE
    $ErrorActionPreference='Stop'
    if($FastpExit-ne 0){throw "fastp failed: $($Group.Name) (exit $FastpExit)"}
}
$Reports=@(Get-ChildItem -LiteralPath $OutRoot -Filter '*.fastp.json')
if($Reports.Count-ne 12){throw "Expected 12 fastp JSON reports; found $($Reports.Count)"}
"PASS`t12/12 fastp reports`t$(Get-Date -Format o)" | Set-Content (Join-Path $RunRoot 'STEP02_FASTP.PASS.txt') -Encoding UTF8
Write-Output 'PASS: GSE193136 raw fastp reports generated for 12/12 samples; no cleaned FASTQ was written'
