# GSE274832 software environment

The accepted run used:

- Windows 10 x64 with PowerShell 5.1
- WSL Linux x86_64
- fastp 1.3.6
- Salmon 1.11.4
- GENCODE v50 / GRCh38.p14
- Python 3.12.13
- pandas 3.0.1
- NumPy 2.3.5
- R 4.5.3
- tximport 1.38.2
- DESeq2 1.50.2
- ggplot2 4.0.2
- pheatmap 1.0.13

The full R `sessionInfo()` from the accepted analysis is stored at:

`06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE274832/R_SESSION_INFO.txt`

The Salmon index build requires an `en_US.UTF-8` locale inside WSL because Salmon 1.11.4 requests that locale during final index serialization.
