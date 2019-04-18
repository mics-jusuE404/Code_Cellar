## Template for ChIP-seq differential analysis using CSAW.
## Run everything default but then use called peaks to aggregate windows:

library(csaw)
library(GenomicAlignments)
library(data.table)
library(edgeR)
library(BSgenome.Mmusculus.UCSC.mm10)

setwd("~/IMTB/Fischer2019/ChIP-seq/")

bam.files <- list.files(path='~/IMTB/Fischer2019/ChIP-seq/bam/', pattern=glob2rx("PU.1*.bam"), full.names = T)
BLACKLIST <- fread("/Volumes/Rumpelkammer/Genomes/mm10/mm10_consensusBL.bed", data.table = F, header = F)
frag.len <- 200
#######################################################################################################

## Estimate fragment length:
max.delay <- 500
x <- correlateReads(bam.files, max.delay, param=param)
plot(0:max.delay, x, type="l", ylab="CCF", xlab="Delay (bp)")

## Define blacklists:
tmp.seqlevels_all <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
tmp.seqlevels_keep <- tmp.seqlevels_all[!grepl("chrM|random|chrU|chrY",tmp.seqlevels_all)]
tmp.blacklist <- GRanges(seqnames = BLACKLIST[,1], ranges = IRanges(start = BLACKLIST[,2]+1, end = BLACKLIST[,3]))
param <- readParam(BPPARAM = MulticoreParam(workers=detectCores()-1), restrict = tmp.seqlevels_keep, discard = tmp.blacklist)

## Read data excluding the blacklisted regions:
data <- windowCounts(bam.files, ext=frag.len, width=10, param=param)

## Keep 1% of the windows assuming (in concordance with peak calling macs2)
## that 1% of the genome is PU.1 bound:
keep <- filterWindows(data, type="proportion")$filter > 0.99
sum(keep)
filtered.data <- data[keep,]

## Normalization:
binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)
filtered.data <- normFactors(binned, se.out=filtered.data)
filtered.data$norm.factors

## Check by MA plot:
adj.counts <- cpm(asDGEList(binned), log=TRUE)
normfacs <- filtered.data$norm.factors
for (i in seq_len(length(bam.files)-1)) {
  cur.x <- adj.counts[,1]
  cur.y <- adj.counts[,1+i]
  smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y,
                xlab="A", ylab="M", main=paste("1 vs", i+1))
  all.dist <- diff(log2(normfacs[c(i+1, 1)]))
  abline(h=all.dist, col="red")
}
                
## edgeR:
y <- asDGEList(filtered.data)
y$samples$group <- c(1,1,0,0)
design <- model.matrix(~1+group, data=y$samples)
rownames(design) <- c("PU.1_Ca_rep1", "PU.1_Ca_rep2", "PU.1_Wt_rep1", "PU.1_Wt_rep2")

## Estimate and plot dispersion:
y <- estimateDisp(y, design)
#par(mfrow=c(1,2))
#o <- order(y$AveLogCPM)
#plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", 
#     lwd=2,ylim=c(0, 1), 
#     xlab=expression("Ave."~Log[2]~"CPM"),
#     ylab=("Biological coefficient of variation"))
#fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

results <- glmQLFTest(fit, contrast=c(0, 1))
rowData(filtered.data) <- cbind(rowData(filtered.data), results$table)

## Merge bins by called peaks:
called_peaks <- fread("bam/PU1_calledCombined_peaks_filtered.narrowPeak", select = c(1,2,3), data.table=F, header = F)
GR <- GRanges(seqnames = called_peaks[,1],
              ranges = IRanges(start = called_peaks[,2],
                               end = called_peaks[,3]))
olap <- findOverlaps(GR, rowRanges(filtered.data))
tabbroad <- combineOverlaps(olap, results$table)

df.final <- cbind(data.frame(GR, data.frame(tabbroad)))
