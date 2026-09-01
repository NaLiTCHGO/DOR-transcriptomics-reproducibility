param(
    [string]$ProjectRoot = '.',
    [string]$BatchRunId = '20260824_BATCH_F3_S6_V04',
    [string]$SourceBatchRunId = '20260824_BATCH_F3_S6_V03'
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$ProjectRoot = (Resolve-Path -LiteralPath $ProjectRoot).Path
$RunsRoot = Join-Path $ProjectRoot '05_analysis_steps\PUB_FIGURE_REBUILD\runs'
$SourceBatchRoot = (Resolve-Path -LiteralPath (Join-Path $RunsRoot $SourceBatchRunId)).Path
$RunRoot = Join-Path $RunsRoot $BatchRunId
$SkillRoot = 'external\research_scientific_figure_suite'
$Rscript = 'C:\Program Files\R\R-4.5.3\bin\Rscript.exe'

if (-not (Test-Path -LiteralPath $SkillRoot -PathType Container)) { throw "Active Skill 3 root is missing: $SkillRoot" }
if (-not (Test-Path -LiteralPath (Join-Path $SkillRoot 'SKILL.md') -PathType Leaf)) { throw 'Active Skill 3 SKILL.md is missing' }
if (-not (Test-Path -LiteralPath $Rscript -PathType Leaf)) { throw "Rscript is missing: $Rscript" }

$FormalScripts = @(
    (Join-Path $ProjectRoot '04_code\R\PUB_FIGBATCH_00_common.R'),
    (Join-Path $ProjectRoot '04_code\R\PUB_FIGBATCH_02_render_main.R'),
    (Join-Path $ProjectRoot '04_code\R\PUB_FIGBATCH_03_render_supp_cohorts.R'),
    (Join-Path $ProjectRoot '04_code\R\PUB_FIGBATCH_04_render_supp_synthesis.R'),
    (Join-Path $ProjectRoot '04_code\R\PUB_FIGBATCH_05_automated_qc.R'),
    (Join-Path $ProjectRoot '04_code\R\PUB_FIGBATCH_06_build_unified_visual_qc.R'),
    (Join-Path $ProjectRoot '04_code\R\PUB_FIGBATCH_07_reconstruct_vst_heatmaps.R'),
    $PSCommandPath
)
foreach ($Script in $FormalScripts) {
    if (-not (Test-Path -LiteralPath $Script -PathType Leaf)) { throw "Formal script is missing: $Script" }
}

$Directories = @('figures\main', 'figures\supplement', 'reconstructed_heatmaps', 'qc', 'manifests', 'logs', 'code_snapshot')
foreach ($RelativePath in $Directories) {
    New-Item -ItemType Directory -Path (Join-Path $RunRoot $RelativePath) -Force | Out-Null
}

if (-not (Test-Path -LiteralPath (Join-Path $RunRoot 'source_data') -PathType Container)) {
    Copy-Item -LiteralPath (Join-Path $SourceBatchRoot 'source_data') -Destination $RunRoot -Recurse
}
foreach ($Script in $FormalScripts) {
    Copy-Item -LiteralPath $Script -Destination (Join-Path $RunRoot 'code_snapshot') -Force
}

# Only Supplementary Figures S1-S3 are rebuilt in this round. Every other
# batch candidate was accepted in the preceding visual-QC round and is copied
# byte-for-byte with its existing warning/provenance records.
$FrozenFiles = @(
    'figures\main\Figure_3_Gene_Reproducibility_and_Heterogeneity_v09_QC.png',
    'figures\main\Figure_4_Pathway_Convergence_v05_QC.png',
    'figures\main\Figure_5_LOCO_Robustness_v08_QC.png',
    'figures\supplement\Figure_S4_Gene_Synthesis_Sensitivity_v05_QC.png',
    'figures\supplement\Figure_S5_Pathway_Context_v07_QC.png',
    'figures\supplement\Figure_S6_LOCO_Context_v04_QC.png',
    'qc\RUNTIME_WARNING_LOG_MAIN.csv',
    'qc\FULL_RUNTIME_WARNING_RECONCILIATION_MAIN.md',
    'qc\RUNTIME_WARNING_LOG_SUPP_SYNTHESIS.csv',
    'qc\FULL_RUNTIME_WARNING_RECONCILIATION_SUPP_SYNTHESIS.md',
    'manifests\RENDERER_PROVENANCE_MAIN.md',
    'manifests\R_SESSION_INFO_MAIN.txt',
    'manifests\RENDERER_PROVENANCE_SUPP_SYNTHESIS.md',
    'manifests\R_SESSION_INFO_SUPP_SYNTHESIS.txt'
)
foreach ($RelativePath in $FrozenFiles) {
    $SourcePath = Join-Path $SourceBatchRoot $RelativePath
    $DestinationPath = Join-Path $RunRoot $RelativePath
    if (-not (Test-Path -LiteralPath $SourcePath -PathType Leaf)) { throw "Frozen source artifact is missing: $SourcePath" }
    Copy-Item -LiteralPath $SourcePath -Destination $DestinationPath -Force
}

$Commands = @"
# User Visual QC Round 4 execution plan

- Batch run: $RunRoot
- Source lock / frozen candidates: $SourceBatchRoot
- Active Skill 3: $SkillRoot
- Scientific claim/result recomputation: NO
- Display-only reconstruction: Figure S1-S3 panel-D Top-30 VST matrices
- New rendering: Figure S1-S3 v08 with native pheatmap column labels
- Frozen without rerender: Figure 3 v09, Figure 4 v05, Figure 5 v08, Figure S4 v05, Figure S5 v07, Figure S6 v04

## Saved-script order

1. PUB_FIGBATCH_07_reconstruct_vst_heatmaps.R
2. PUB_FIGBATCH_03_render_supp_cohorts.R
3. PUB_FIGBATCH_05_automated_qc.R
4. PUB_FIGBATCH_06_build_unified_visual_qc.R

All formal scripts were copied to code_snapshot/ before execution. The reconstruction is accepted only if normalized counts match the locked tables and clustered sample order matches the locked panel.
"@
Set-Content -LiteralPath (Join-Path $RunRoot 'RUN_COMMANDS.md') -Value $Commands -Encoding utf8

$FrozenCard = @"
# Frozen-candidate carry-forward card

- Source batch: $SourceBatchRoot
- Destination batch: $RunRoot
- Figure 3-5 and Supplementary Figures S4-S6 were copied without rerendering.
- Reason: the user found no remaining visual defect in those batch figures in round 4.
- Scientific claim/result recomputation: NO
"@
Set-Content -LiteralPath (Join-Path $RunRoot 'FROZEN_CANDIDATE_CARRY_FORWARD.md') -Value $FrozenCard -Encoding utf8

function Write-Status {
    param([string]$State, [string]$Step, [string]$Message)
    $Body = @"
# User Visual QC Round 4 status

- State: $State
- Current/completed step: $Step
- Message: $Message
- Batch run: $RunRoot
- Updated: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss zzz')
"@
    Set-Content -LiteralPath (Join-Path $RunRoot 'RUN_STATUS.md') -Value $Body -Encoding utf8
}

function Invoke-RStage {
    param([string]$Stage, [string]$ScriptName, [string[]]$ExtraArguments)
    $ScriptPath = Join-Path $ProjectRoot "04_code\R\$ScriptName"
    $LogPath = Join-Path $RunRoot "logs\$Stage.log"
    Write-Status -State 'RUNNING' -Step $Stage -Message "Executing saved script $ScriptPath"
    $PreviousErrorActionPreference = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    & $Rscript --vanilla $ScriptPath "--project-root=$ProjectRoot" "--run-root=$RunRoot" "--skill-root=$SkillRoot" @ExtraArguments 2>&1 |
        Tee-Object -FilePath $LogPath
    $ExitCode = $LASTEXITCODE
    $ErrorActionPreference = $PreviousErrorActionPreference
    if ($ExitCode -ne 0) {
        Write-Status -State 'FAIL' -Step $Stage -Message "Exit code $ExitCode; see $LogPath"
        throw "$Stage failed with exit code $ExitCode"
    }
}

$PreviousLcAll = $env:LC_ALL
$PreviousLang = $env:LANG
$env:LC_ALL = 'C'
$env:LANG = 'C'
try {
    Invoke-RStage -Stage 'STEP01_RECONSTRUCT_S1_S3_VST_DISPLAYS' -ScriptName 'PUB_FIGBATCH_07_reconstruct_vst_heatmaps.R' -ExtraArguments @()
    Invoke-RStage -Stage 'STEP02_RENDER_S1_S3' -ScriptName 'PUB_FIGBATCH_03_render_supp_cohorts.R' -ExtraArguments @()
    Invoke-RStage -Stage 'STEP03_AUTOMATED_QC' -ScriptName 'PUB_FIGBATCH_05_automated_qc.R' -ExtraArguments @()
    Invoke-RStage -Stage 'STEP04_UNIFIED_VISUAL_QC' -ScriptName 'PUB_FIGBATCH_06_build_unified_visual_qc.R' -ExtraArguments @()
    Write-Status -State 'PASS' -Step 'STEP04_UNIFIED_VISUAL_QC' -Message 'Round-4 candidates rendered and automated QC passed; assistant actual-render review is next.'
}
finally {
    if ($null -eq $PreviousLcAll) { Remove-Item Env:LC_ALL -ErrorAction SilentlyContinue } else { $env:LC_ALL = $PreviousLcAll }
    if ($null -eq $PreviousLang) { Remove-Item Env:LANG -ErrorAction SilentlyContinue } else { $env:LANG = $PreviousLang }
}

Write-Output "PASS: user visual-QC round 4 candidate batch completed: $RunRoot"
