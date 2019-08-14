## Rscript to calculate TMM size factors from a count matrix read from stdin,
## which is then used in the bigwig script to normalize bigwig counts.

## turn off scientific notation to avoid rounding errors:
options(scipen=99999)

## check for presence of required packages:
packageS <- c("data.table", "edgeR")
if (length(grep("FALSE", (packageS %in% rownames(installed.packages())))) > 0){
  stop("Package(s): ", packageS[which( packageS %in% rownames(installed.packages()) == "FALSE")], " are not installed!")
}

suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(data.table))

## read data:
raw.counts <- fread('cat /dev/stdin', skip = 1, header = T, data.table = F)

if (ncol(raw.counts) == 7) {
  write.table(data.frame(c("ONESAMPLE")), sep="\t", col.names = F, row.names = T,
              quote = F, file = "./sizeFactors.txt")
  stop("Only one sample present, so skipping calculation!")              
}

## expecting featureCounts format remove non-count columns
raw.counts <- raw.counts[,7:ncol(raw.counts)]

## normFactors:
tmp.NF <- calcNormFactors(object = raw.counts, method = c("TMM"))

## raw library size:
tmp.LS <- colSums(raw.counts)

## effective normalization factor => each count in the bigwig will be divided by this
if (file.exists("./TMMfactors.txt")) {
  Do.Append <- TRUE
} else Do.Append <- FALSE

suppressWarnings(
  write.table(x = data.frame(Sample    = names(tmp.LS),
                             TMMfactor = tmp.NF * ( tmp.LS / 1000000)),
              file = paste0("./TMMfactors.txt"),
              append = Do.Append, sep = "\t",
              col.names = !Do.Append, row.names = FALSE, quote = FALSE)
)
