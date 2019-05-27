# Content:

## DNAseq_lowlevel_vx.sh
This script performs all relevant preprocessing of DNA-seq data (currently with settings for ChIP-seq and ATAC-seq) starting with `fastqc` control, adapter trimming, alignment, filtering, peak calling, FRiP calculation, library complexity assessment, browser track creation, summary statistics/multiqc. THere are a few options that can be set in the header such as number of parallel jobs and presets for ATACseq and ChIPseq and the input format and layout (uBam or Fastq, paired- or single-end runs, but not mixed). Expects input files to be in same directory as the script itself.

## ATAC-seq
Folder contains some scripts to analyze ATAC-seq data such as a wrapper to run `chromVAR`.

## ChIP-seq
Same as ATAC-seq but for, guess what, ChIP-seq data.

## HiC
Draft for a wrapper for the HiCUPs pipeline.

## RNA-seq
Wrappers for quantifications with `salmon` or `hisat2`. `hisat2` script also contains some quality assessment via RSeqC. 
Still, `salmon` is the preferred quantification pipeline as it respects multimappers and implements GC bias correction with the method implemented in the BioC package `alpine` plus is directly compatible with `tximport` to create depth/length offsets for use in `edgeR`.

## scRNA-seq
Draft script for quantification of scRNA-seq reads per cell via `alevin`. Currently more orless untested and preliminary.

## Differential_Analysis
Folder contains wrapper scripts to run `csaw` (ATAC-seq/ChIP-seq), `edgeR` (RNA-seq) given `tximport` input, `tximport` to import transcript abundance estimated from `salmon` and some scripts to plot volcanos and MA-plots in a more or less automated/reproducible fashion.
