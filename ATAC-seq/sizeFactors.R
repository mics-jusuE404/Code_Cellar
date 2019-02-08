#!/home/a/a_toen03/anaconda3_things/anaconda3/envs/R_env/bin/Rscript

## Script to calculate DESeq2 size factors for ATAC-seq data based on a count matrix
## over the entire genome with 500bp windows, taking the top 100k windows based on rowMeans.
## File comes from stdin via the main bash script ATACseq_lowlevel(...).sh

suppressMessages(library(DESeq2))
suppressMessages(library(data.table))
 
## read data:
counts <- fread('cat /dev/stdin', skip = 1, header = T, data.table = F)

## get rowMeans
counts <- data.frame(counts, rowMean = rowMeans(counts[,7:ncol(counts)]))

## sort by rowMeans
counts <- counts[with(counts, order(-rowMean)),]

## top 100k regions based on rowMeans:
counts_subset <- head(counts[,7:(ncol(counts)-1)], n=100000)

## size factors:
SizeFactors <- estimateSizeFactorsForMatrix(counts = counts_subset)

## write:
write.table(data.frame(SizeFactors), sep="\t", col.names = F, row.names = T,
            quote = F, file = "./sizeFactors.txt")
