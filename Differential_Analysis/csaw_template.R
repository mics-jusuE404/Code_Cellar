## Template for ATAC-seq differential analysis using CSAW:

library(csaw)
library(GenomicAlignments)
library(data.table)
library(edgeR)

setwd("~/IMTB/Fischer2019/ChIP-seq/")

bam.files <- list.files(path='~/IMTB/Fischer2019/ChIP-seq/bam/', pattern=glob2rx("PU.1*.bam"), full.names = T)
BLACKLIST <- fread("/Volumes/Rumpelkammer/Genomes/mm10/mm10_consensusBL.bed", data.table = F, header = F)
frag.len <- 200

peaks <- fread("peaks/pu1_merged.bed", data.table = F)
peaks.gr <- GRanges(seqnames = peaks[,1], ranges=IRanges(start=peaks[,2]+1, end=peaks[,3]))

#######################################################################################################

## Estimate fragment length:
#max.delay <- 500
#x <- correlateReads(bam.files, max.delay, param=param)
#plot(0:max.delay, x, type="l", ylab="CCF", xlab="Delay (bp)")

## Define blacklists:
tmp.seqlevels_all <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
tmp.seqlevels_keep <- tmp.seqlevels_all[!grepl("chrM|random|chrU|chrY",tmp.seqlevels_all)]
tmp.blacklist <- GRanges(seqnames = BLACKLIST[,1], ranges = IRanges(start = BLACKLIST[,2]+1, end = BLACKLIST[,3]))
param <- readParam(BPPARAM = MulticoreParam(workers=detectCores()-1), restrict = tmp.seqlevels_keep, discard = tmp.blacklist)

## Read data excluding the blacklisted regions:
data <- regionCounts(bam.files, ext=frag.len, param=param, regions = peaks.gr)

## Normalization:
binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)
data <- normFactors(data, se.out=data)

plotMA_smooth <- function(counts, main){
  R=counts[,1]
  G=counts[,2]
  M=log2(R/G)
  A=0.5*log2(R*G)
  smoothScatter(A, M, xlab="A", ylab="M", main = main); abline(h=0)
}

plotMA_smooth(counts = assay(data)[,c(1,2)], "raw")
plotMA_smooth(counts = calculateCPM(data, use.norm.factors = T, log = F)[,c(1,2)], "norm")

## edgeR:
y <- asDGEList(data)
y$samples$group <- c(1,1,0,0)
design <- model.matrix(~1+group, data=y$samples)
rownames(design) <- c("PU.1_Ca_rep1", "PU.1_Ca_rep2", "PU.1_Wt_rep1", "PU.1_Wt_rep2")

## Estimate and plot dispersion:
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
results <- glmQLFTest(fit, contrast=c(0, 1))
out <- topTags(results, n=Inf, adjust.method="BH", sort.by="none")

## Combine ranges and results (careful, only works of topTags is not sorted!)
out.final <- data.frame( data.frame(rowRanges(data))[,1:3], out$table)
colnames(out.final)[1] <- "chr"


