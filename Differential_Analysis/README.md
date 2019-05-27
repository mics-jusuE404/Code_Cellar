# Differential_Analysis

## csaw:
Wrapper to perform differential analysis of DNA-seq count data such as ATAC-seq or ChIP-seq using `csaw` and `edgeR`.
Currently, it accepts a list of peak summit positions that are then extended to a user-defined size, followed by read counting directly from the BAM files. It then normalizes data using TMM using either the counts in peaks or counts in large bins of 10kb (the latter as suggested in the `csaw` manual). Default is peak counts as it usually (for ATAC-seq, based on examination of MA-plots) better scales data under the assumption that most regions are not differential. Wrapper also outputs MA-plots for either all sample combinations or averaged over replicate groups to check normalization efficiency. 

## edgeR:
Differential RNA-seq analysis using the `edgeR` framework downstream of `salmon` and `tximport`. Expects `txi` from `tximport` plus design matrix and contrasts as inputs. Outputs the standard `topTags` from `edgeR` plus the offset-corrected CPMs (offsets that account for depth, composition bias and txlength from `tximport`).

## plotMA_custom:
A wrapper to create MA-plots for pairwise comparison of count data based on average counts and log2FC.

## tximport:
Short template script to summarize `salmon` quantifications to the gene level.

## volcano:
Simple wrapper to create colored volcano plots based on p-values and logFCs.
