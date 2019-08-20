## Wrapper script for InTAD analysis, taking several args from the command line.
## This version is still tailored to our current project and not generic as it expects
## very specific input from the --rdata flag.

##############################################################################################################################

## Helper to stop quietly without the ugly stop message from stop(), cudos to Stibu
## https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
SQ <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

##############################################################################################################################

## Check if required packages are present:
message("[Info]: Checking / loading packages")
packageS <- c("InTAD", "data.table", "parallel", "mclust", "optparse")

if (length(grep("FALSE", (packageS %in% rownames(installed.packages())))) > 0){
  message(
    paste("[Error]:", 
          "Package(s): [",
          paste(packageS[which( packageS %in% rownames(installed.packages()) == "FALSE")], collapse = " "),
          "] are not installed!"))
  SQ()
}

for (i in packageS){
  suppressPackageStartupMessages(require(i, character.only = T))
}

##############################################################################################################################

## required helper function:
source("makeUniqueCombinations.R")
source("findCutoff_mclustBIC.R")
source("findCorrelation_custom.R")

list.sources <- c("makeUniqueCombinations", "findCutoff_mclustBIC", "findCorrelation_custom")

if ( sum(grepl(paste(paste0("^", list.sources, "$"), collapse = "|"), ls())) 
     !=
     length(list.sources)){
  message("[Error]: Some of the required functions from source could not be loaded, check that!")
  SQ()
}

##############################################################################################################################

## Arg parsing:
option_list = list(
  
  make_option(c("--study"), action="store", default=NA, type="character",
              help="The study name to filter the count matrix for (character)"),
  
  make_option(c("--cores"), action="store", default=floor(parallel::detectCores()/2), type="integer",
              help="Number of cores to work on one shuffling. Default is half of available cores (integer)"),
  
  make_option(c("--nshuff"), action="store", default=1, type="integer",
              help="Number of shufflings to perform. Default is 1 (integer)"),
  
  make_option(c("--nshuffParallel"), action="store", default=1, type="integer",
              help="Number of parallel shuffling processes. Default is 1 (integer)"),
  
  make_option(c("--rdata"), action="store", default=NA, type="character",
              help="Path to the .Rdata file to load containing input data"),
  
  make_option(c("--outdir"), action="store", default="./InTAD_results", type="character",
              help="The directory to store results in, default is ./InTAD_results")  
  
); opt = parse_args(OptionParser(option_list=option_list))

arg.mis <- c()
Study          = opt$study; if(is.na(Study)) arg.mis <- c(arg.mis, ("--study"))
Cores          = opt$cores
Shuffs         = opt$nshuff
nshuffParallel = opt$nshuffParallel
Rdata          = opt$rdata; if(is.na(Rdata)) arg.mis <- c(arg.mis, ("--rdata"))
Outdir         = opt$outdir

if (!dir.exists(Outdir)) dir.create(Outdir, recursive = T)

## Stop if some are NA empty:
if (length(arg.mis) > 0) {
  message("[Error]: Those options are missing without defaults:")
  for (i in arg.mis){ message(paste0("         ", i))}
  SQ()
  message("         Run the script with the --help flag to get usage information")
}

if (Cores > parallel::detectCores()) {
  Cores <- parallel::detectCores()
  message(paste("[Info]: Setting --cores to", parallel::detectCores(), "as this is the maximum on this node"))
}

if ( (Cores * nshuffParallel) > parallel::detectCores()){
  nshuffParallel <- floor(parallel::detectCores()/Cores)
  message(paste("[Info]: Setting --nshuffParallel to", 
                nshuffParallel,
                "as original --nshuffParallel would've exceeded number of available workers"))
}

if (nshuffParallel > 2) message("[Warning]: --nshuffParallel is > 2, should only be done on the largesmp nodes!")

##############################################################################################################################

## Load the environment with the preprocessed input data
message("[Info]: Loading environment")
load(paste0(Rdata))
if (length(ls()) == 0) {
  message("[Error]: Environment is empty, check if --rdata object contains correct input data")
  SQ()
}

##############################################################################################################################

## set explicit seed 
rm(.Random.seed)
set.seed(2019)

##############################################################################################################################

## Prepare data from --rdata (specific for this project, should be made more generic in the future)
counts.atac = atac.fpm_clusters[,grep(Study, colnames(atac.fpm_clusters))]
ranges.atac = atac.ranges_clusters.gr
counts.rna = rna.fpkm[,grep(Study, colnames(rna.fpkm))]
ranges.rna = rna.ranges
TADs = tads.gr

##############################################################################################################################

## Make table with all possible sample combinations to draw from (=Shuffs), then randomize:
message("[Info]: Shuffling samples to get all possible combinations")
permut.combos <- makeUniqueCombinations(list.rna  = colnames(counts.rna),
                                        list.atac = colnames(counts.atac))
permut.combos <- sample(permut.combos, length(permut.combos))
message(paste("[Info]: Obtained", length(permut.combos), "possible combinations"))

## save to disk for reproducibility:
write.table(data.frame(Number=seq(1,length(permut.combos)),
                       Name=permut.combos),
            file = paste0(Outdir,
                          "/InTAD_shuffles_",
                          Study,
                          ".tsv"),
            col.names = T,
            row.names = F,
            sep="\t",
            quote = F)

if (length(permut.combos) < Shuffs){
  message(
    paste("[Info]: Setting --nshuffs to", length(permut.combos), "as this is the maximum available number"))
  Shuffs <- length(permut.combos)
}

##############################################################################################################################

## Find a gene expression cutoff:
message("[Info]: Estimating expression cutoff")
cutVal <- round(findCutoff_mclustBIC(count.matrix = counts.rna, do.log2 = TRUE, log2.prior = 1, plotExprDistr = FALSE), 
                digits=3)

message(paste("[Info]: Filtering RNA-seq data for at least one group-mean above (based on log2+1)", cutVal))

tmp.unique.celltype <- unique(sapply(strsplit(colnames(counts.rna), split="_rep"), function(m) m[1]))

keep.rna <-
  unique(
    unlist(
      lapply(tmp.unique.celltype, function(U){
        tmp.which <- grep(paste0("^",U, "_rep"), colnames(counts.rna))
        return(as.integer(which(rowMeans(counts.rna[,tmp.which]) >= cutVal)))
      }
      )
    )
  )

tmp.number <- nrow(counts.rna)
counts.rna <- counts.rna[keep.rna,]
ranges.rna <- ranges.rna[keep.rna,]

message("[Info]: From initially ", tmp.number, " transcripts retained ", nrow(counts.rna), " after filtering ")

## Note, as we use TSS positions so transcript counts rather than gene counts we have a lot of tests.
## It might be worth to eventually feed only one tss per gene to qvalue(), maybe the one that has the best correlation
## to avoid high mt burden, have to think about that.

##############################################################################################################################

## mclapply function for parallelization around the main workhorse functions,
## use 2 parallel on 1.5Tb nore or 4 on the 3Tb node, not more!
message("[Info]: Starting InTAD analysis")
message("")
be.quiet <- mclapply(seq(1, Shuffs), mc.cores = nshuffParallel, function(SHUFF){
  
  ############################################################################################################################
  
  ## Select the samples for the current shuffling:
  tmp.split <- strsplit(strsplit(permut.combos[SHUFF], split="##|@@")[[1]], split = "__")
  where.rna <- c()
  where.atac <- c()
  for (i in 1:length(tmp.split)){
    tmp.rna.extract  <- tmp.split[[i]][1]
    tmp.atac.extract <- tmp.split[[i]][2]
    where.rna  <- c(where.rna,  grep(paste0("^", tmp.rna.extract, "$"), colnames(counts.rna)))
    where.atac <- c(where.atac, grep(paste0("^", tmp.atac.extract, "$"), colnames(counts.atac)))
  }
  
  ## select the chosen samples:
  tmp.counts.rna   <- counts.rna [,where.rna]
  tmp.counts.atac  <- counts.atac[,where.atac]
  
  ## rename to ensure paired samples have matching names:  
  colnames(tmp.counts.rna)  <- paste0("sample", seq(1,length(where.rna)))
  colnames(tmp.counts.atac) <- paste0("sample", seq(1,length(where.atac)))
  
  ############################################################################################################################
  
  ## Make InTADSig object:
  InSigTad <- newSigInTAD(signalData = tmp.counts.atac, signalRegions = ranges.atac,
                          countsData = tmp.counts.rna, geneRegions = ranges.rna, performLog = TRUE, 
                          logExprsOffset = 1, ncores = Cores)
  
  ## Combine in TAD:
  InSigTad <- combineInTAD(object = InSigTad, tadGR = tads.gr, selMaxTadOvlp = FALSE)
  
  ## results w/o qvalues, this will be done later as we might want to filter out some things before
  corRes   <- findCorrelation_custom(object = InSigTad, method = "pearson", adj.pval = FALSE, plot.proportions = FALSE)
  
  ## write results to disk:
  write.table(corRes, file = paste0(Outdir,
                                    "/InTAD_results_",
                                    Study,
                                    "_shuffle",
                                    SHUFF,
                                    ".tsv"),
              col.names = T, row.names = F, sep="\t", quote = F)
  
}
) ## end mclapply
