## Template for differential analysis with csaw using called peaks as template for read counting and FDR control:

library(csaw)
library(GenomicAlignments)
library(data.table)

###################################################################################################################

## Script assumes that 
setwd("~/IMTB/Fischer2019/ChIP-seq/")

###################################################################################################################
###################################################################################################################

## Wrapper to run csaw in peakbased mode:
run_csaw_peakBased <- function(NAME="", SUMMITS, WIDTH=200, 
                               BAMS, FRAGLEN="", PAIRED=F, 
                               NORM="largebins", GROUP, CONTRASTS, 
                               CORES=c(detectCores()-2)){
  
  #################################################################################################################
  message("")
  message("[STARTED]: ", NAME)
  
  if ( (NAME == "") || (CONTRASTS == "") || (GROUP == "") ){
    stop("Please enter a {NAME} for the analysis -- Exiting")
    break
  }
  
  if (PAIRED == T) {
    stop("Paired-end mode not yet implemented -- Exiting")
    break
  }

 
  #################################################################################################################
  ## Take summits and extend to desired window size
  if ( class(SUMMITS) != "GRanges") {
    tmp.gr <- GRanges(seqnames = SUMMITS[,1],
                      ranges = IRanges(start = SUMMITS[,2]+1, end = SUMMITS[,3]))
  }
  
  ## Resize and collapse in case of overlaps:
  message("Extending peak summits to ", WIDTH, "bp")
  peaks_resized.gr <- GenomicRanges::reduce(GenomicRanges::resize(tmp.gr, width = WIDTH, fix = "center"))
  
  #################################################################################################################
  ## Fragment length:
  #if (PAIRED == T) {
    
    #PSE="both"
    #out <- getPESizes(pe.bam)
  #}
    
  if (PAIRED == F) PSE="none"
  
  PARAM=readParam(BPPARAM = MulticoreParam(workers = CORES), 
                  discard = makeGRangesFromDataFrame(df = BLACKLIST, 
                                                     seqnames.field = "V1",
                                                     start.field = "V2",
                                                     end.field = "V3"), 
                  pe = PSE)
  
  if (FRAGLEN > 0) message("Using an average fragment length of ", FRAGLEN, "bp")
    
  if (PAIRED == F & FRAGLEN == ""){
    
    message("Calculating average fragment length")
    max.delay <- 500
    x <- correlateReads(BAMS, max.delay, param=PARAM)
    FRAGLEN <- maximizeCcf(x)
    
    ## Save xcorrplot:
    if (FRAGLEN > 0){
      plot.new()
      assign( paste(NAME, "_plot_xcorr", sep=""), 
              recordPlot(plot(0:max.delay, x, type="l", ylab="CCF", xlab="Delay (bp)")),
              envir = .GlobalEnv)
    }; dev.off()
    
  } 
  
  #################################################################################################################
  ## Read BAM:
  message("Counting reads per peak")
  data <- regionCounts(BAMS, 
                       ext = FRAGLEN, 
                       param = PARAM,
                       regions = peaks_resized.gr)
  
  #################################################################################################################
  ## Normalize
  if (NORM == "largebins"){
    message("Calculating normalization factors")
    binned <- windowCounts(BAMS, bin=TRUE, width=10000, param=PARAM)
    data <- normFactors(binned, se.out=data)
  }
  
  if (NORM == "peaks") {
    data <- normFactors(data, se.out=data)
  }
  
  ## Save counts and put in colnames:
  strReverse <- function(x){
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
  }
  file.names <- as.character(strReverse(sapply(strReverse(BAMS), function(x) strsplit(x, split = "\\/")[[1]][1])))
  
  RAWcounts <- data.frame(assay(data)); colnames(RAWcounts) <- file.names
  CPMcounts <- data.frame(calculateCPM(data, log = F)); colnames(CPMcounts) <- file.names
  
  ## Save counts:
  assign( paste(NAME, "_data_RAWcounts", sep=""),
          RAWcounts, envir = .GlobalEnv)
  assign( paste(NAME, "_data_CPMcounts", sep=""),
          CPMcounts, envir = .GlobalEnv)
  
  ## MA plots for average over the conditions:
  #assign( paste(NAME, "_plot_xcorr", sep=""), 
  #        #MAPLOT code
  #        envir = .GlobalEnv)
  
  #################################################################################################################
  ## edgeR part:
  message("Running edgeR")
  
  ## set up DGElist:
  y <- asDGEList(data)
  y$samples$group <- GROUP
  
  ## design matrix:  
  design <- model.matrix(~ 0 + group, data=y$samples)
  colnames(design) <- gsub("group", "", colnames(design))
  assign( paste(NAME, "_data_design", sep=""),
          design, envir = .GlobalEnv)
  
  ## estimate dispersion:
  message("Estimating dispersion")
  y <- estimateDisp(y, design)
  assign( paste(NAME, "_data_estimateDisp", sep=""),
          y, envir = .GlobalEnv)
  
  ## fit the model:
  message("Fitting GLM")
  fit <- glmQLFit(y, design, robust=TRUE)
  assign( paste(NAME, "_data_fit", sep=""),
          fit, envir = .GlobalEnv)
  
  ## test all contrasts:
  message("Testing all contrasts")
  for (i in CONTRASTS){
    
    ## current contrast:
    current.results <- glmQLFTest(fit, contrast = makeContrasts(contrasts = i, levels=design))
    assign( paste(NAME, "_data_resultsRAW_", gsub("-", "_",i), sep=""),
            current.results, envir = .GlobalEnv)
    
    ## fdr:
    current.out <- topTags(current.results, n=Inf, adjust.method="BH", sort.by="none")
    current.out <- data.frame( data.frame(rowRanges(data))[,1:3], current.out$table)
    assign( paste(NAME, "_data_resultsFDR_", gsub("-", "_",i), sep=""),
            current.out, envir = .GlobalEnv)

  }
  
  #################################################################################################################

}

###################################################################################################################
###################################################################################################################


## Example:
BLACKLIST <- fread("/Volumes/Rumpelkammer/Genomes/mm10/mm10_consensusBL.bed", data.table = F, header = F)
PU1_Summits <- fread("peaks/PU.1_combinedCall_summits.bed", data.table = F)  
PU1_Bams <- list.files(path='~/IMTB/Fischer2019/ChIP-seq/bam/', pattern=glob2rx("PU.1*.bam"), full.names = T)

run_csaw_peakBased(NAME = "PU1", SUMMITS = PU1_Summits, 
                   WIDTH = 200,  BAMS = PU1_Bams, 
                   PAIRED = F,   GROUP = c("Ca", "Ca", "Wt", "Wt"), 
                   CONTRASTS = c("Ca-Wt"), FRAGLEN = "", 
                   NORM = "largebins", CORES=14)
