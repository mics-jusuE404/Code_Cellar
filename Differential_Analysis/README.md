# Differential_Analysis

## edgeR_DNAseq_v1.R
Script for differential analysis of ChIP- and ATAC-eq experiments starting from a count matrix in GRanges format with raw data.
Will perform count normalization, plotting of MA-plots (to check norm.), PCA, and saving CPMs to disk.
If `Design=NULL` it will stop after doing the above. If a design is provided, it will proceed with differential analysis,
testing the provided contrasts, optionally against a certain FC. Will then plot MAs with signif (5%FDR) regions colored in red.
Eventually will save the analysis as Rdata.
If one has a count matrix from `featureCounts` use:

`counts.gr <- makeGRangesFromDataFrame(counts[,c(2,3,4,c(relevant_sample_columns)], keep.extra.columns = T)`

By default the script will trim suffixes `_dedup.bam` as this is our standard convention for readily-filtered BAM files.
Can be changed in the function options. Parameters re explained on top of the script.

## edgeR_RNAseq_v1.R
Pretty much the same as the DNAseq script but for, guess what, RNA-seq. Starts from salmon quantifications and then performs all the stuff as above.

## smoothScatter_TopTags.R
Takes a topTags object from edgeR, produces an MA-plot using the `smoothScatter()` function and colors significant regions (FDR < threshold) in red.
Example:
![MAcolor](https://i.ibb.co/ccDytB7/Screen-Shot-2019-05-29-at-19-31-33.png)



