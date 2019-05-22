# Bash_Bioinformatics

## ATAC-seq
Folder contains the low-level script to process ATAC-seq data including quality control, adapter trimming, alignment, peak calling, FRiP calculation, bigwig creation, preseq library complexity and multiqc summary. Also contains a wrapper for running `chromVAR` to infer accessability changes linked to transcription factor motif occurrence.

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
