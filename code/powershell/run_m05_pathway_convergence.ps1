param(
    [string]$ProjectRoot = "",
    [string]$RunId = "20260814_M05_B1_MSIGDB2026_1",
    [string]$RscriptExe = "C:\Program Files\R\R-4.5.3\bin\Rscript.exe",
    [switch]$Force
)

$ErrorActionPreference = "Stop"

if ([string]::IsNullOrWhiteSpace($ProjectRoot)) {
    $ProjectRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot "..\.."))
} else {
    $ProjectRoot = [System.IO.Path]::GetFullPath($ProjectRoot)
}

$runDir = Join-Path $ProjectRoot "05_analysis_steps\M05_PATHWAY_CONVERGENCE\runs\$RunId"
$logDir = Join-Path $runDir "logs"
$passMarker = Join-Path $runDir "M05_PASS_WITH_LIMITATION.txt"
$completeMarker = Join-Path $runDir "RUN_COMPLETE.txt"
$analysisScript = Join-Path $ProjectRoot "04_code\R\M05_01_run_ranked_pathway_convergence.R"
$auditScript = Join-Path $ProjectRoot "04_code\python\M05_02_audit_pathway_convergence.py"

New-Item -ItemType Directory -Path $runDir -Force | Out-Null
New-Item -ItemType Directory -Path $logDir -Force | Out-Null

if ((Test-Path -LiteralPath $passMarker) -and (Test-Path -LiteralPath $completeMarker) -and -not $Force) {
    Write-Output "M05 already completed: $runDir"
    exit 0
}

foreach ($required in @($RscriptExe, $analysisScript, $auditScript)) {
    if (-not (Test-Path -LiteralPath $required)) { throw "Required M05 executable/script missing: $required" }
}

$python = (Get-Command python -ErrorAction Stop).Source
$started = Get-Date

Write-Output "[M05 01/02] Cohort-first ranked fgsea and pathway synthesis"
$savedErrorAction = $ErrorActionPreference
$ErrorActionPreference = "Continue"
& $RscriptExe --vanilla $analysisScript --project-root $ProjectRoot --run-dir $runDir 1> (Join-Path $logDir "01_fgsea.stdout.log") 2> (Join-Path $logDir "01_fgsea.stderr.log")
$fgseaExit = $LASTEXITCODE
$ErrorActionPreference = $savedErrorAction
if ($fgseaExit -ne 0) { throw "M05 fgsea failed with exit code $fgseaExit" }

Write-Output "[M05 02/02] Independent output audit"
& $python $auditScript --run-dir $runDir 1> (Join-Path $logDir "02_audit.stdout.log") 2> (Join-Path $logDir "02_audit.stderr.log")
if ($LASTEXITCODE -ne 0) { throw "M05 output audit failed with exit code $LASTEXITCODE" }

if (-not (Test-Path -LiteralPath $passMarker)) { throw "M05 audit exited successfully but PASS marker is missing" }

$finished = Get-Date
$completeText = @(
    "module=M05_PATHWAY_CONVERGENCE"
    "run_id=$RunId"
    "status=PASS_WITH_LIMITATION"
    "started=$($started.ToString('yyyy-MM-dd HH:mm:ss zzz'))"
    "finished=$($finished.ToString('yyyy-MM-dd HH:mm:ss zzz'))"
    "elapsed_minutes=$([math]::Round(($finished - $started).TotalMinutes, 2))"
    "run_dir=$runDir"
) -join [Environment]::NewLine
[System.IO.File]::WriteAllText($completeMarker, $completeText, [System.Text.UTF8Encoding]::new($false))

Write-Output "M05 completed with PASS_WITH_LIMITATION: $runDir"
