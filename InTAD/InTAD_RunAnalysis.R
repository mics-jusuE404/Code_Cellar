## Customized (and in parts ugly) InTAD wrapper script intended to be run on the 1.5TB or 3TB RAM node
## It assumes all required files (scripts, Rdata) to be in the same directory.

##############################################################################################################################

## Check for packages:
packageS <- c("InTAD", "data.table", "parallel", "mclust")
if (length(grep("FALSE", (packageS %in% rownames(installed.packages())))) > 0){
  stop("Package(s): ", packageS[which( packageS %in% rownames(installed.packages()) == "FALSE")], " are not installed!")
}

## load packages:
suppressMessages(require(parallel))
suppressMessages(require(data.table))
suppressMessages(require(InTAD))
suppressMessages(require(mclust))

##############################################################################################################################

## Load the environment with the preprocessed input data
load("190807_InTAD_PreparedData.Rdata")
if (length(ls()) == 0) stop("Environment appears empty, load() probably not successful!")

##############################################################################################################################

## set explicit seed 
rm(.Random.seed)
set.seed(2019)

##############################################################################################################################

## load the customized functions that were modified from the original InTAD source from BioC
source("makeUniqueCombinations.R")
source("findCorrelation_custom.R")
source("findCutoff_mclustBIC.R")

## Some dummy data for testing purposes
Do.Test <- FALSE
if (Do.Test){
  Study = "Rasmussen"
  counts.atac = atac.fpm_clusters[seq(1,10),grep(Study, colnames(atac.fpm_clusters))]
  ranges.atac = atac.ranges_clusters.gr[seq(1,10)]
  counts.rna = rna.fpkm[,grep(Study, colnames(rna.fpkm))]
  ranges.rna = rna.ranges
  TADs = tads.gr
  permut.cycles = 1 
  ExprThresh = 2
  Studyname = Study
  do.log2 = TRUE 
  log2.prior = 1
}

##############################################################################################################################

## Wrapper script for the customized InTAD analysis:
## => processing signal pairs in chunks of 1000, perform shuffling of RNA/ATAC samples
InTAD_wrapper <- function(counts.atac,                 ## fpkm-norm. ATAC-seq counts
                          ranges.atac,                 ## associated coordiantes
                          counts.rna,                  ## same for RNA-seq as above
                          ranges.rna, 
                          TADs,                        ## TAD coordinates
                          permut.cycles = 1,           ## number of permutations, default no shuffling
                          chunk.size = 10000,           ## size of chunks to split cor.test
                          ExprThresh = NULL,           ## an expression threshold for transcript/gene expression
                          ## we calculate this externally using mclust
                          FilterByGroupAverage = TRUE, ## keep only genes when the celltype average is above ExprThresh
                          ## FALSE would be to have at least one sample above this threshold
                          Studyname,                   ## the study name of the samples to be considered
                          ## necessary as the count matrices include different studies
                          do.log2 = TRUE,              ## do everything on log2-transformed counts
                          log2.prior = 1){             ## prior count for log2 transformation to avoid log2(0)
  
  if(is.na(Studyname)) stop("Something wrong with studyname read from command line!")
  
  ##############################################################################################################################
  ##
  ## Some background info:
  ## Typically when using mouse data you do not have matched ATAC-seq and RNA-seq from the exact same donor,
  ## but from different littermates. Therefore when matching samples one has no guarantee the result is actually
  ## 100% representative. One therefore should perform multiple correlation analysis, each time shuffling the
  ## RNA-seq and ATAC-seq samples of the same cell type to reduce the impact of combination-specific significances
  ## The rationale is that a correlation (averaged or combined by e.g. Fisher's Method based on many shufflings)
  ## is more trustworthy than a single correlation result. Therefore perform multiple correlation analysis.
  ##
  ##############################################################################################################################
  
  ## For this, take the input data and create a list with all possible combinations of the samples,
  ## respecting (=not mixing up) cell types. So far this is only implemented for two different cell types:
  message("[Info]: Shuffling samples to get all possible combinations")
  permut.combos <- makeUniqueCombinations(list.rna  = colnames(counts.rna),
                                          list.atac = colnames(counts.atac))
  message(paste("[Info]: Obtained", length(permut.combos), "possible combinations"))
  
  ## randomize them:
  permut.combos <- sample(permut.combos, length(permut.combos))
  
  if (permut.cycles > length(permut.combos)) {
    message(paste("[Info]: Setting number of shufflings to", length(permut.combos),
                  "as requested number", paste0("(",permut.cycles,")"),
                  "was higher than available combinations"))
    permut.cycles <- length(permut.combos)
  }
  
  ##############################################################################################################################
  
  ## Estimate cutoff for non-expressed genes:
  if (do.log2) counts.rna <- log2(counts.rna + log2.prior)
  
  message("[Info]: Estimating RNA-seq expression cutoff")
  rna.cutoff <- findCutoff_mclustBIC(count.matrix = counts.rna, 
                                     do.log2 = FALSE,
                                     plotExprDistr = FALSE)
  
  rna.cutoff <- round(rna.cutoff, digits = 3)
  
  if (!is.numeric(rna.cutoff)) stop("Something went wrong with RNA-seq cutoff calculation!")
  
  ##########
  ## Case 1: FilterByGroupAverage == TRUE, so filter for genes where the average group (= cell type)
  ##         expression is above this threshold for at least one group of cells.
  ##         This is more conservative than the default in InTAD (which is at least one sample above thresh.),
  ##         but in this specific case (littermate wildtype RNA-seq) where variation is typically
  ##         limited probably a good idea as it reduces the candidate TSS by factor > 2,
  ##         reducing memory consumption and multiple testing burden, while probably not really
  ##         filtering out biolgically-meaningful genes anyway.
  ##
  ##
  if (FilterByGroupAverage){
    
    tmp.unique.celltype <- unique(sapply(strsplit(colnames(counts.rna), split="_rep"), function(m) m[1]))
    
    keep.rna <-
      unique(
        unlist(
          lapply(tmp.unique.celltype, function(U){
            tmp.which <- grep(paste0("^",U, "_rep"), colnames(counts.rna))
            return(as.integer(which(rowMeans(counts.rna[,tmp.which]) >= rna.cutoff)))
          }
          )
        )
      )
  }
  
  ##########
  ## Case 2: Filter for TSS where at least one sample in any grup is above threshold
  if (!FilterByGroupAverage){
    keep.rna <- which(as.numeric(apply(counts.rna, 1, function(x) max(x))) >= rna.cutoff)
  }
  
  ## Filter matrix and coordiantes for selected TSS
  message(paste("[Info]: Filtering RNA-seq data with cutoff =", rna.cutoff))
  counts.rna <- counts.rna[keep.rna,]
  ranges.rna <- ranges.rna[keep.rna,]
  
  ############################################################################################################################
  ########## Now run the actual InTAD analysis:
  
  trash.can <- lapply(seq(1,permut.cycles), function(Q){
    
    ## split the string that indicates how to combine the ATAC/RNA pairs,
    ## see the makeUniqueCombinations.R for details:
    tmp.split <- strsplit(strsplit(permut.combos[Q], split="##|@@")[[1]], split = "__")
    
    ## find the positions of the elements that are to be considered in this run in the count matrices:
    where.rna <- c()
    where.atac <- c()
    
    ## keep only the selected samples in the respective combination:  
    for (i in 1:length(tmp.split)){
      tmp.rna.extract  <- tmp.split[[i]][1]
      tmp.atac.extract <- tmp.split[[i]][2]
      where.rna  <- c(where.rna,  grep(paste0("^", tmp.rna.extract, "$"), colnames(counts.rna)))
      where.atac <- c(where.atac, grep(paste0("^", tmp.atac.extract, "$"), colnames(counts.atac)))
    }
    
    ## select the chosen samples:
    tmp.rna   <- counts.rna [,where.rna]
    tmp.atac  <- counts.atac[,where.atac]
    
    ## rename to ensure paired samples have matching names:  
    colnames(tmp.rna) <- paste0("sample", seq(1,length(where.rna)))
    colnames(tmp.atac) <- paste0("sample", seq(1,length(where.atac)))
    
    ##########################################################################################################################
    ## Start of actual InTAD part, the above was data sanitation and preparation:
    ## Combine everything into the InTAD wrapper
    message("[Info]: Making InTadSig object")
    tmp.intad <- newSigInTAD(signalData = tmp.atac, 
                             signalRegions = ranges.atac, 
                             countsData    = tmp.rna, 
                             geneRegions   = ranges.rna,
                             performLog = TRUE)
    
    message("[Info]: Combining genes and signals per TAD")                  
    tmp.intad <- combineInTAD(object = tmp.intad, 
                              tadGR = TADs, 
                              selMaxTadOvlp = FALSE)
    
    ## Stop if no connections were found
    if (nrow(tmp.intad@signalConnections[[1]]) == 0){
      stop("[Error]: Did not find any valid candidate gene/enhancer pairs!")
    }
    
    ## As it turned out that the actual analysis is extremely memory-consuming for large
    ## numbers of signal/TSS combinations we will by default split it into chunks.
    ## For this we use the modified findCorrelation.R function
    ## It will split the cor analysis into chunks of (default) 5000 signals at a time using a lapply-based function
    ## and spill the result to disk folder "./InTAD_raw",
    ## using 72 workers for the correlation analysis itself 
    ## If you use more workers use the fat servers (1.5 or 3TB RAM) instead of the 192GB ones
    trashy.can <- findCorrelation_custom(object = tmp.intad, 
                                         Studyname = Studyname,
                                         method = "pearson",
                                         chunk.size = chunk.size,
                                         current.cycle = Q,           ## Q is the current variable from the lapply above
                                         total.cycle = permut.cycles  
    )
    }
    )
}

## Accept args for Study name from command line:
CLargs = commandArgs(trailingOnly=TRUE)
Study  = CLargs[1]

## Now the actual command to be executed.
## Script itself is to be called via Rscript thisscript.R <studyname> 2> stderr.log
InTAD_wrapper(counts.atac = atac.fpm_clusters[,grep(Study, colnames(atac.fpm_clusters))], 
              ranges.atac = atac.ranges_clusters.gr, 
              counts.rna = rna.fpkm[,grep(Study, colnames(rna.fpkm))], 
              ranges.rna = rna.ranges, 
              TADs = tads.gr, 
              permut.cycles = CLargs[2], ## number of shufflings to perform
              ExprThresh = 2,            ## leave hardcoded for this analysis
              Studyname = Study,     ##study name to filter count matrix for
              do.log2 = TRUE, 
              log2.prior = 1,
              chunk.size = 25000) ## careful, only got that high on the big nodes with 1.5 or 3Tb RAM!

##############################################################################################################################
