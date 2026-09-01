param(
    [string]$ProjectRoot = "",
    [string]$RunId = "20260814_M06_B1_PATHWAY_LOCO",
    [switch]$Force
)

$ErrorActionPreference = "Stop"
if ([string]::IsNullOrWhiteSpace($ProjectRoot)) {
    $ProjectRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot "..\.."))
} else {
    $ProjectRoot = [System.IO.Path]::GetFullPath($ProjectRoot)
}

$runDir = Join-Path $ProjectRoot "05_analysis_steps\M06_LEAVE_ONE_COHORT_OUT\runs\$RunId"
$logDir = Join-Path $runDir "logs"
$passMarker = Join-Path $runDir "M06_PASS_WITH_LIMITATION.txt"
$completeMarker = Join-Path $runDir "RUN_COMPLETE.txt"
$analysisScript = Join-Path $ProjectRoot "04_code\python\M06_01_run_pathway_loco.py"
$auditScript = Join-Path $ProjectRoot "04_code\python\M06_02_audit_pathway_loco.py"
New-Item -ItemType Directory -Path $runDir -Force | Out-Null
New-Item -ItemType Directory -Path $logDir -Force | Out-Null

if ((Test-Path -LiteralPath $passMarker) -and (Test-Path -LiteralPath $completeMarker) -and -not $Force) {
    Write-Output "M06 already completed: $runDir"
    exit 0
}
foreach ($required in @($analysisScript, $auditScript)) {
    if (-not (Test-Path -LiteralPath $required)) { throw "Required M06 script missing: $required" }
}
$python = (Get-Command python -ErrorAction Stop).Source
$started = Get-Date
Write-Output "[M06 01/02] Rotation-wise pathway selection and held-out replication"
& $python $analysisScript --project-root $ProjectRoot --run-dir $runDir 1> (Join-Path $logDir "01_loco.stdout.log") 2> (Join-Path $logDir "01_loco.stderr.log")
if ($LASTEXITCODE -ne 0) { throw "M06 pathway LOCO failed with exit code $LASTEXITCODE" }
Write-Output "[M06 02/02] Independent no-leakage and output audit"
& $python $auditScript --run-dir $runDir 1> (Join-Path $logDir "02_audit.stdout.log") 2> (Join-Path $logDir "02_audit.stderr.log")
if ($LASTEXITCODE -ne 0) { throw "M06 audit failed with exit code $LASTEXITCODE" }
if (-not (Test-Path -LiteralPath $passMarker)) { throw "M06 audit passed but marker is missing" }
$finished = Get-Date
$completeText = @(
    "module=M06_LEAVE_ONE_COHORT_OUT"
    "run_id=$RunId"
    "status=PASS_WITH_LIMITATION"
    "started=$($started.ToString('yyyy-MM-dd HH:mm:ss zzz'))"
    "finished=$($finished.ToString('yyyy-MM-dd HH:mm:ss zzz'))"
    "elapsed_minutes=$([math]::Round(($finished - $started).TotalMinutes, 2))"
    "run_dir=$runDir"
) -join [Environment]::NewLine
[System.IO.File]::WriteAllText($completeMarker, $completeText, [System.Text.UTF8Encoding]::new($false))
Write-Output "M06 completed with PASS_WITH_LIMITATION: $runDir"
