param(
    [string]$ProjectRoot = "",
    [string]$RunId = "20260814_M04_B1_THREE_CORE_RNASEQ",
    [switch]$Force
)

$ErrorActionPreference = "Stop"

if ([string]::IsNullOrWhiteSpace($ProjectRoot)) {
    $ProjectRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot "..\.."))
} else {
    $ProjectRoot = [System.IO.Path]::GetFullPath($ProjectRoot)
}

$runDir = Join-Path $ProjectRoot "05_analysis_steps\M04_REPRO_HETEROGENEITY\runs\$RunId"
$logDir = Join-Path $runDir "logs"
$passMarker = Join-Path $runDir "M04_PASS_WITH_LIMITATION.txt"
$completeMarker = Join-Path $runDir "RUN_COMPLETE.txt"
$analysisScript = Join-Path $ProjectRoot "04_code\python\M04_01_run_cross_cohort_synthesis.py"
$auditScript = Join-Path $ProjectRoot "04_code\python\M04_02_audit_cross_cohort_results.py"

New-Item -ItemType Directory -Path $runDir -Force | Out-Null
New-Item -ItemType Directory -Path $logDir -Force | Out-Null

if ((Test-Path -LiteralPath $passMarker) -and (Test-Path -LiteralPath $completeMarker) -and -not $Force) {
    Write-Output "M04 already completed: $runDir"
    exit 0
}

foreach ($required in @($analysisScript, $auditScript)) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "Required formal script missing: $required"
    }
}

$python = (Get-Command python -ErrorAction Stop).Source
$started = Get-Date

Write-Output "[M04 01/02] Cross-cohort synthesis"
& $python $analysisScript --project-root $ProjectRoot --run-dir $runDir 1> (Join-Path $logDir "01_synthesis.stdout.log") 2> (Join-Path $logDir "01_synthesis.stderr.log")
if ($LASTEXITCODE -ne 0) {
    throw "M04 synthesis failed with exit code $LASTEXITCODE"
}

Write-Output "[M04 02/02] Independent output audit"
& $python $auditScript --run-dir $runDir 1> (Join-Path $logDir "02_audit.stdout.log") 2> (Join-Path $logDir "02_audit.stderr.log")
if ($LASTEXITCODE -ne 0) {
    throw "M04 output audit failed with exit code $LASTEXITCODE"
}

if (-not (Test-Path -LiteralPath $passMarker)) {
    throw "M04 audit exited successfully but PASS marker is missing"
}

$finished = Get-Date
$completeText = @(
    "module=M04_REPRO_HETEROGENEITY"
    "run_id=$RunId"
    "status=PASS_WITH_LIMITATION"
    "started=$($started.ToString('yyyy-MM-dd HH:mm:ss zzz'))"
    "finished=$($finished.ToString('yyyy-MM-dd HH:mm:ss zzz'))"
    "elapsed_minutes=$([math]::Round(($finished - $started).TotalMinutes, 2))"
    "run_dir=$runDir"
) -join [Environment]::NewLine
[System.IO.File]::WriteAllText($completeMarker, $completeText, [System.Text.UTF8Encoding]::new($false))

Write-Output "M04 completed with PASS_WITH_LIMITATION: $runDir"
