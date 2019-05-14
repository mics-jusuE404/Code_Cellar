## Template for differential analysis with csaw using called peaks as template for read counting and FDR control:

run_csaw_peakBased <- function(NAME="", SUMMITS, WIDTH=200, 
                               BAMS, FRAGLEN="", PAIRED=F, 
                               NORM="largebins", GROUP, CONTRASTS, 
                               CORES=c(detectCores()-2), UNREP=F,
                               RERUN=F){
  
  library(csaw)
  library(edgeR)
  library(GenomicAlignments)
  library(data.table)
  
  #################################################################################################################
  message("")
  message("[STARTED]: ", NAME)
  
  ## message if unreplicated design:
  if (UNREP == T) message("[INFO] Running on unreplicated data!")
  
  ## message if no name given in function
  if ( (NAME == "") || (CONTRASTS == "") || (GROUP == "") ){
    stop("Please enter a {NAME} for the analysis -- Exiting")
    break
  }
  
  ## stop if paired end specified as not impemented yet
  if (PAIRED == T) {
    stop("Paired-end mode not yet implemented -- Exiting")
    break
  }
  
  ## RERUN=F means run everything as usual, T would be to skip all counting part and do edgeR right away.
  ## can be useful if MA plots show that normalization was not good and one simply wants to repeat edgeR with 
  ## the existing counts but a different normalization
  if (RERUN == F) {
    assign( paste(NAME, "rerun_check", sep=""), c(1), envir = .GlobalEnv)
  }
  
  if (RERUN == T) {
    if ( ! exists( paste(NAME, "rerun_check", sep="") ) ) {
      stop("RERUN was set to TRUE but function has not been run before")
    } else {
      message("Skipping upstream part and starting with edgeR part right away")
    }
  }
      
  #################################################################################################################
  if (RERUN == F) {
    
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
  
  ## PARAM, discarding blacklisted regions
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
    
    assign(paste(NAME, "_fraglength", sep=""), FRAGLEN, envir = .GlobalEnv)
    
    ## Save xcorrplot:
    if (FRAGLEN > 0){
      plot(0:max.delay, x, type="l", ylab="CCF", xlab="Delay (bp)")
      assign( paste(NAME, "_plot_xcorr", sep=""), 
              recordPlot(),
              envir = .GlobalEnv)
    }
    
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
  message("Calculating normalization factors")
  binned <- windowCounts(BAMS, bin=TRUE, width=10000, param=PARAM)
  
  ## Save counts and put in colnames:
  strReverse <- function(x){
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
  }
  file.names <- as.character(strReverse(sapply(strReverse(BAMS), function(x) strsplit(x, split = "\\/")[[1]][1])))
  
  RAWcounts <- data.frame(assay(data))
  colnames(RAWcounts) <- file.names
  
  CPMcountsPeaks <- data.frame( calculateCPM(object = normFactors(data, se.out=data), use.norm.factors = T, log = F) )
  colnames(CPMcountsPeaks) <- file.names
  
  CPMcountsBins  <- data.frame( calculateCPM(object = normFactors(binned, se.out=data), use.norm.factors = T, log = F) )
  colnames(CPMcountsBins) <- file.names
  
  if (NORM == "peaks") {
    data <- normFactors(data, se.out=data)
  }
  
  if (NORM == "largebins") {
    data <- normFactors(binned, se.out=data)
  }
  
  ## Save counts:
  assign( paste(NAME, "_data_countsRAW", sep=""),
          RAWcounts, envir = .GlobalEnv)
  
  assign( paste(NAME, "_data_countsCPM_peaks", sep=""),
          CPMcountsPeaks, envir = .GlobalEnv)
  
  assign( paste(NAME, "_data_countsCPM_largebins", sep=""),
          CPMcountsBins, envir = .GlobalEnv)
  
  } ## END OF RERUN
  
  ## Extra part for RERUN=T:
  if (RERUN == T){
    if (NORM == "peaks"){
      message("RERUN set to TRUE -- using peak-derived normalization factors")
      data <- normFactors(data, se.out=data)
    }
    if (NORM == "largebins"){
      message("RERUN set to TRUE -- using largebin-derived normalization factors")
      data <- normFactors(binned, se.out=data)
    }
  }
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
  
  ## Model fitting:
  if (UNREP == F){
    
    ## estimate dispersion:
    message("Estimating dispersion")
    y <- estimateDisp(y, design)
    assign( paste(NAME, "_data_estimateDisp", sep=""),
            y, envir = .GlobalEnv)
    
    message("Fitting GLM")
    fit <- glmQLFit(y, design, robust=TRUE)
    assign( paste(NAME, "_data_fit", sep=""),
            fit, envir = .GlobalEnv)
    
  }
  
  ## If unreplicated:
  if (UNREP == T){
    message("Skipping dispersion estimates and use 0.1 as proxy")
    message("Fitting GLM on UNREPLICATED data with disp = 0.1")
    fit <- glmFit(y, design, dispersion=0.1)
    assign( paste(NAME, "_data_Unrepfit", sep=""),
            fit, envir = .GlobalEnv)
  }
  
  #################################################################################################################
  ## test all contrasts:
  message("Testing all contrasts")
  for (i in CONTRASTS){
    
    if (UNREP == F){
      
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
    
    if (UNREP == T){
      
      ## LRT tests:
      current.results <- glmLRT(fit, contrast = makeContrasts(contrasts = i, levels=design))
      assign( paste(NAME, "_data_resultsUnrepRAW_", gsub("-", "_",i), sep=""),
              current.results, envir = .GlobalEnv)
      
      ## fdr:
      current.out <- topTags(current.results, n=Inf, adjust.method="BH", sort.by="none")
      current.out <- data.frame( data.frame(rowRanges(data))[,1:3], current.out$table)
      assign( paste(NAME, "_data_resultsUnrepFDR_", gsub("-", "_",i), sep=""),
              current.out, envir = .GlobalEnv) 
    }

  }
  
  #################################################################################################################
  
  message("[ENDED]: ", NAME)
  message("")
  
}
