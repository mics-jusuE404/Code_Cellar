# Differential_Analysis

## csaw:
Wrapper to perform differential analysis of DNA-seq count data such as ATAC-seq or ChIP-seq using `csaw` and `edgeR`.
Currently, it accepts a list of peak summit positions that are then extended to a user-defined size, followed by read counting directly from the BAM files. It then normalizes data using TMM using either the counts in peaks or counts in large bins of 10kb (the latter as suggested in the `csaw` manual). Default is peak counts as it usually (for ATAC-seq, based on examination of MA-plots) better scales data under the assumption that most regions are not differential. Wrapper also outputs MA-plots for either all sample combinations or averaged over replicate groups to check normalization efficiency. 

## salmon2edgeR:
Wrapper to aggregate `salmon` quantifications to the gene level followed by count normalization with `edgeR` using the length-offsets from `tximport`. By default script uses TX2Gene file on my local machine for mouse gencode annotations, and also by default will filter out genes that code for smallRNA as I think these RNAs are poorly captured in standard RNA-seq, therefore sources of false-positives and will inflate multiple testing burden. Wrapper will save the CPM-normalized counts as `"./Lists/*_CPM.tsv"` and also produce some MDS and MA-plots to get a basic idea on how the data look like. The `*.DGElist` that the wrapper produces in `R` can then be used for DEG analysis in `edgeR` starting with the `filterByExpr()` function.

## smoothScatter_TopTags.R
Pretty much the same as `plotMA_custom` but takes a topTags object and colors significant regions. 
Example:
![MAcolor](https://i.ibb.co/ccDytB7/Screen-Shot-2019-05-29-at-19-31-33.png)

## volcano:
Simple wrapper to create colored volcano plots based on p-values and logFCs.


