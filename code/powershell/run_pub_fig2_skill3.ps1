param(
    [string]$ProjectRoot = '.',
    [string]$RunId = '20260824_FIG2_SKILL3_V04'
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path -LiteralPath $ProjectRoot).Path
$skillRoot = 'external\research_scientific_figure_suite'
$rscript = 'C:\Program Files\R\R-4.5.3\bin\Rscript.exe'
$runDir = Join-Path $project "05_analysis_steps\PUB_FIGURE_REBUILD\runs\$RunId"
$formalScripts = @(
    (Join-Path $project '04_code\R\PUB_FIG2_01_prepare_source_data.R'),
    (Join-Path $project '04_code\R\PUB_FIG2_02_render_candidate.R'),
    (Join-Path $project '04_code\R\PUB_FIG2_03_automated_qc.R'),
    (Join-Path $project '04_code\R\PUB_FIG2_04_hotspot_and_font_qc.R'),
    (Join-Path $project '04_code\powershell\run_pub_fig2_skill3.ps1')
)

if (-not (Test-Path -LiteralPath $skillRoot)) { throw "Active Skill 3 root not found: $skillRoot" }
if (-not (Test-Path -LiteralPath (Join-Path $skillRoot 'SKILL.md'))) { throw 'Active Skill 3 SKILL.md is missing' }
if (-not (Test-Path -LiteralPath $rscript)) { throw "Rscript not found: $rscript" }
foreach ($script in $formalScripts) {
    if (-not (Test-Path -LiteralPath $script)) { throw "Formal script missing: $script" }
}

$dirs = @('source_data', 'source_data\upstream_refs', 'figures', 'legends', 'qc', 'manifests', 'logs', 'failed_attempts', 'code_snapshot')
foreach ($relative in $dirs) {
    $target = Join-Path $runDir $relative
    if (-not (Test-Path -LiteralPath $target)) { New-Item -ItemType Directory -Path $target | Out-Null }
}

foreach ($script in $formalScripts) {
    Copy-Item -LiteralPath $script -Destination (Join-Path $runDir ('code_snapshot\' + [IO.Path]::GetFileName($script))) -Force
}

$commands = @"
# Figure 2 Skill 3 run commands

- Run ID: $RunId
- Active Skill: $skillRoot
- R: $rscript --vanilla
- Scientific boundary: locked M02/M03 tables only; no PCA/DESeq2 recomputation.

## Executed order

1. PUB_FIG2_01_prepare_source_data.R --project-root PROJECT --run-dir RUN
2. PUB_FIG2_02_render_candidate.R --project-root PROJECT --run-dir RUN --skill-root SKILL
3. PUB_FIG2_03_automated_qc.R --run-dir RUN
4. PUB_FIG2_04_hotspot_and_font_qc.R --run-dir RUN

All five formal scripts were copied into code_snapshot/ before the first R step.
"@
Set-Content -LiteralPath (Join-Path $runDir 'RUN_COMMANDS.md') -Value $commands -Encoding utf8

function Write-RunStatus {
    param([string]$State, [string]$Step, [string]$Message)
    $body = @"
# Figure 2 run status

- Run ID: $RunId
- State: $State
- Current/completed step: $Step
- Message: $Message
- Run folder: $runDir
- Updated: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss zzz')
"@
    Set-Content -LiteralPath (Join-Path $runDir 'RUN_STATUS.md') -Value $body -Encoding utf8
}

function Invoke-RStep {
    param([string]$Name, [string]$Script, [string[]]$Arguments)
    $log = Join-Path $runDir ("logs\" + $Name + '.log')
    Write-RunStatus -State 'RUNNING' -Step $Name -Message "Executing saved script $Script"
    & $rscript --vanilla $Script @Arguments 2>&1 | Tee-Object -FilePath $log
    $exitCode = $LASTEXITCODE
    if ($exitCode -ne 0) {
        $failure = @"
# Failed attempt: $Name

- Exit code: $exitCode
- Script: $Script
- Log: $log
- Status: FAIL_NOT_AUTHORITATIVE
"@
        Set-Content -LiteralPath (Join-Path $runDir ("failed_attempts\" + $Name + '.md')) -Value $failure -Encoding utf8
        Write-RunStatus -State 'FAIL' -Step $Name -Message "Exit code $exitCode; see failed_attempts and log."
        throw "$Name failed with exit code $exitCode"
    }
}

$previousLcAll = $env:LC_ALL
$previousLang = $env:LANG
$env:LC_ALL = 'C'
$env:LANG = 'C'
try {
    Invoke-RStep -Name 'STEP01_SOURCE_PREPARATION' -Script $formalScripts[0] -Arguments @('--project-root', $project, '--run-dir', $runDir)
    Invoke-RStep -Name 'STEP02_RENDER_V07' -Script $formalScripts[1] -Arguments @('--project-root', $project, '--run-dir', $runDir, '--skill-root', $skillRoot)
    Invoke-RStep -Name 'STEP03_AUTOMATED_QC' -Script $formalScripts[2] -Arguments @('--run-dir', $runDir)
    Invoke-RStep -Name 'STEP04_HOTSPOT_AND_FONT_QC' -Script $formalScripts[3] -Arguments @('--run-dir', $runDir)
    Write-RunStatus -State 'PASS' -Step 'STEP04_HOTSPOT_AND_FONT_QC' -Message 'Source lock, v07 Other-key collision repair, automated QC, hotspot crops and font preflight passed; strict human actual-render review is next.'
}
finally {
    if ($null -eq $previousLcAll) { Remove-Item Env:LC_ALL -ErrorAction SilentlyContinue } else { $env:LC_ALL = $previousLcAll }
    if ($null -eq $previousLang) { Remove-Item Env:LANG -ErrorAction SilentlyContinue } else { $env:LANG = $previousLang }
}

Write-Output "PASS: Figure 2 Skill 3 v07 candidate run completed: $runDir"
