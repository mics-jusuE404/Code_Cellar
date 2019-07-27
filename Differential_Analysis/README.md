# Differential_Analysis

## edgeR_DNAseq_v1.R
Script for differential analysis of ChIP- and ATAC-eq experiments starting from a count matrix in GRanges format with raw data.
Will perform count normalization, plotting of MA-plots (to check norm.), PCA, and saving CPMs to disk.
If `Design=NULL` it will stop after doing the above. If a design is provided, it will proceed with differential analysis,
testing the provided contrasts, optionally against a certain FC. Will then plot MAs with signif (5%FDR) regions colored in red.
Eventually will save the analysis as Rdata.
If one has a count matrix from `featureCounts` use:

`counts.gr <- makeGRangesFromDataFrame(counts[,c(2,3,4,c(relevant_sample_columns)], keep.extra.columns = T)` 

## salmon2edgeR:
Wrapper to aggregate `salmon` quantifications to the gene level followed by count normalization with `edgeR` using the length-offsets from `tximport`. By default script uses TX2Gene file on my local machine for mouse gencode annotations, and also by default will filter out genes that code for smallRNA as I think these RNAs are poorly captured in standard RNA-seq, therefore sources of false-positives and will inflate multiple testing burden. Wrapper will save the CPM-normalized counts as `"./Lists/*_CPM.tsv"` and also produce some MDS and MA-plots to get a basic idea on how the data look like. The `*.DGElist` that the wrapper produces in `R` can then be used for DEG analysis in `edgeR` starting with the `filterByExpr()` function. Can be followed by `run_edgeR.R` for DGE.

## edgeR_RNAseq_V1.R
RNA-seq workflow with `edgeR` assuming input from `salmon2edgeR`.

## smoothScatter_TopTags.R
Pretty much the same as `plotMA_custom` but takes a topTags object and colors significant regions. 
Example:
![MAcolor](https://i.ibb.co/ccDytB7/Screen-Shot-2019-05-29-at-19-31-33.png)



