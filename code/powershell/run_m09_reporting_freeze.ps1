param(
    [Parameter(Mandatory=$true)][string]$ProjectRoot,
    [string]$NodeExe = 'node',
    [string]$PythonExe = 'python',
    [switch]$Force
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path -LiteralPath $ProjectRoot).Path
$run = Join-Path $project '05_analysis_steps\M09_REPORTING_FREEZE\runs\20260815_M09_B1_FINAL_FREEZE'
$runtime = Join-Path $project '00_project_control\artifact_runtime_m09'
$output = Join-Path $project 'outputs\DOR_M09_20260815'
$delivery = Join-Path $project '10_delivery_packages\DOR_Final_Scientific_Package_20260815.zip'
$workbook = Join-Path $output 'DOR_Project_Final_Overview_20260815.xlsx'
$previews = Join-Path $run 'previews'
$results = Join-Path $run 'results'
$logs = Join-Path $run 'logs'

if ((Test-Path -LiteralPath (Join-Path $run 'RUN_COMPLETE.txt')) -and -not $Force) {
    throw 'M09 already completed. Use -Force only for an intentional full rerun.'
}
New-Item -ItemType Directory -Force -Path $runtime,$output,$previews,$results,$logs,(Split-Path $delivery -Parent) | Out-Null
Copy-Item -LiteralPath (Join-Path $project '04_code\javascript\M09_01_build_final_overview.mjs') -Destination (Join-Path $runtime 'M09_01_build_final_overview.mjs') -Force

Write-Host '[M09 01/02] Build and render final overview workbook'
& $NodeExe (Join-Path $runtime 'M09_01_build_final_overview.mjs') --projectRoot $project --output $workbook --previewDir $previews 1> (Join-Path $logs '01_workbook.stdout.log') 2> (Join-Path $logs '01_workbook.stderr.log')
if ($LASTEXITCODE -ne 0) { throw "Workbook build failed with exit code $LASTEXITCODE" }

Write-Host '[M09 02/02] Build Markdown/ZIP and audit release'
& $PythonExe (Join-Path $project '04_code\python\M09_02_build_release_and_audit.py') --project-root $project --workbook $workbook --output-dir $results --zip $delivery 1> (Join-Path $logs '02_release.stdout.log') 2> (Join-Path $logs '02_release.stderr.log')
if ($LASTEXITCODE -ne 0) { throw "Release build/audit failed with exit code $LASTEXITCODE" }

Copy-Item -LiteralPath $workbook -Destination $results -Force
"module=M09_REPORTING_FREEZE`nrun_id=20260815_M09_B1_FINAL_FREEZE`nstatus=PASS`ncompleted_at=$(Get-Date -Format o)`nworkbook=$workbook`nzip=$delivery" | Set-Content -LiteralPath (Join-Path $run 'RUN_COMPLETE.txt') -Encoding utf8
Write-Host "M09 completed: $run"

