param(
  [string]$ProjectRoot=$env:DOR_PROJECT_ROOT,
  [string]$ProxyUrl=$env:DOR_HTTP_PROXY
)
$ErrorActionPreference='Stop'
if(-not $ProjectRoot){$ProjectRoot=(Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path}
function To-Wsl([string]$Path){'/mnt/'+$Path.Substring(0,1).ToLower()+'/'+$Path.Substring(3).Replace('\','/')}
$RefRoot=Join-Path $ProjectRoot '03_data\references\GENCODE_v50_GRCh38p14'
$LogRoot=Join-Path $ProjectRoot '05_analysis_steps\M02_COHORT_QC\runs\20260813_M02_B1_GSE274832\logs\reference'
New-Item -ItemType Directory -Force -Path $RefRoot,$LogRoot|Out-Null
$Files=@(
  @{Name='gencode.v50.transcripts.fa.gz';Url='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.transcripts.fa.gz'},
  @{Name='GRCh38.p14.genome.fa.gz';Url='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh38.p14.genome.fa.gz'},
  @{Name='gencode.v50.annotation.gtf.gz';Url='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz'}
)
$Downloads=@()
foreach($F in $Files){
  $Dest=Join-Path $RefRoot $F.Name
  if(-not(Test-Path $Dest)){
    $Temp="$Dest.download"
    if(Test-Path $Temp){
      $WT=To-Wsl $Temp
      $ErrorActionPreference='Continue';wsl.exe --exec gzip -t $WT 2>$null;$GzipExit=$LASTEXITCODE;$ErrorActionPreference='Stop'
      if($GzipExit-eq 0){Move-Item $Temp $Dest -Force;continue}
      Remove-Item $Temp -Force
    }
    $Out=Join-Path $LogRoot ("{0}.stdout.log" -f $F.Name);$Err=Join-Path $LogRoot ("{0}.stderr.log" -f $F.Name)
    $Args=@('--silent','--show-error','--location','--fail','--retry','100','--retry-all-errors','--retry-delay','5','--connect-timeout','60','--speed-limit','1024','--speed-time','180','--output',$Temp,$F.Url)
    if($ProxyUrl){$Args=@('--proxy',$ProxyUrl)+$Args}
    $P=Start-Process curl.exe -WindowStyle Hidden -ArgumentList $Args -RedirectStandardOutput $Out -RedirectStandardError $Err -PassThru
    $Downloads+=[pscustomobject]@{Name=$F.Name;Dest=$Dest;Temp=$Temp;Process=$P}
  }
}
foreach($D in $Downloads){
  $D.Process.WaitForExit()
  if($D.Process.ExitCode-ne 0){
    $WT=To-Wsl $D.Temp
    $ErrorActionPreference='Continue';wsl.exe --exec gzip -t $WT 2>$null;$GzipExit=$LASTEXITCODE;$ErrorActionPreference='Stop'
    if($GzipExit-ne 0){throw "Reference download failed $($D.Name)"}
  }
}
foreach($D in $Downloads){$W=To-Wsl $D.Temp;wsl.exe --exec gzip -t $W;if($LASTEXITCODE-ne 0){throw "Reference gzip failed $($D.Name)"};Move-Item -LiteralPath $D.Temp -Destination $D.Dest -Force}
foreach($F in $Files){$Dest=Join-Path $RefRoot $F.Name;$W=To-Wsl $Dest;wsl.exe --exec gzip -t $W;if($LASTEXITCODE-ne 0){throw "Reference gzip failed $($F.Name)"}}
$Genome=Join-Path $RefRoot 'GRCh38.p14.genome.fa.gz';$Tx=Join-Path $RefRoot 'gencode.v50.transcripts.fa.gz';$Gentrome=Join-Path $RefRoot 'gencode.v50.gentrome.fa.gz';$Decoys=Join-Path $RefRoot 'decoys.txt';$Index=Join-Path $RefRoot 'salmon_index_v1.11.4'
if(-not(Test-Path $Gentrome)){$GentromeTemp="$Gentrome.building";$WTx=To-Wsl $Tx;$WGenome=To-Wsl $Genome;$WGentromeTemp=To-Wsl $GentromeTemp;wsl.exe --exec bash -lc "cat <(gzip -cd '$WTx') <(gzip -cd '$WGenome') | gzip -c > '$WGentromeTemp'";if($LASTEXITCODE-ne 0){throw 'Gentrome creation failed'};wsl.exe --exec gzip -t $WGentromeTemp;if($LASTEXITCODE-ne 0){throw 'Gentrome gzip validation failed'};Move-Item $GentromeTemp $Gentrome -Force}
if(-not(Test-Path $Decoys)){$WGenome=To-Wsl $Genome;$WDecoys=To-Wsl $Decoys;wsl.exe --exec bash -lc "gzip -cd '$WGenome' | grep '^>' | cut -d ' ' -f 1 | sed 's/>//' > '$WDecoys'"}
$Salmon=To-Wsl (Join-Path $ProjectRoot '00_project_control\tools\linux_x86_64\salmon-1.11.4\bin\salmon');$WG=To-Wsl $Gentrome;$WD=To-Wsl $Decoys;$WI=To-Wsl $Index
if(-not(Test-Path (Join-Path $Index 'versionInfo.json'))){
  $IndexTemp="$Index.building"
  if(Test-Path $IndexTemp){$Failed="$Index.failed_"+(Get-Date -Format yyyyMMdd_HHmmss);Move-Item $IndexTemp $Failed}
  $LinuxBuild='/root/dor_salmon_index_v1.11.4_'+(Get-Date -Format yyyyMMdd_HHmmss)+'.building'
  wsl.exe --exec mkdir -p $LinuxBuild
  if($LASTEXITCODE-ne 0){throw 'Linux-local index staging failed'}
  # Salmon 1.11.4 links a formatting facet that explicitly requests
  # en_US.UTF-8 during final index serialization.  A smoke test is retained in
  # tests/salmon_locale_smoke_test.fa; this locale must exist in WSL.
  wsl.exe --exec env LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 $Salmon index -t $WG -d $WD -i $LinuxBuild -p 8 --gencode --keepDuplicates
  if($LASTEXITCODE-ne 0){throw 'Salmon index failed'}
  wsl.exe --exec test -f "$LinuxBuild/versionInfo.json"
  if($LASTEXITCODE-ne 0){throw 'Salmon index completion marker missing'}
  New-Item -ItemType Directory -Force -Path $IndexTemp|Out-Null
  $WIndexTemp=To-Wsl $IndexTemp
  wsl.exe --exec bash -lc "cp -a '$LinuxBuild'/\. '$WIndexTemp/'"
  if($LASTEXITCODE-ne 0){throw 'Copying completed Salmon index into project failed'}
  if(-not(Test-Path (Join-Path $IndexTemp 'versionInfo.json'))){throw 'Copied Salmon index completion marker missing'}
  if(Test-Path $Index){Remove-Item $Index -Recurse -Force}
  Move-Item $IndexTemp $Index
}
Write-Output 'PASS: GENCODE v50 / GRCh38.p14 decoy-aware Salmon index ready'
