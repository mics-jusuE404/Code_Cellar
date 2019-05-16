## Template for differential analysis with csaw using called peaks as template for read counting and FDR control:
## Wrapper to run csaw in peakbased mode:

run_csaw_peakBased <- function(NAME="", SUMMITS, WIDTH=200, 
                               BAMS, FRAGLEN="", PAIRED=F, 
                               NORM="peaks", GROUP, CONTRASTS, 
                               CORES=c(detectCores()-1), 
                               MAall = T, WD = "~/"){
  
  packageS <- c("csaw", "edgeR", "GenomicAlignments", "data.table")
  if (length(grep("FALSE", (packageS %in% rownames(installed.packages())))) > 0){
    stop("Package(s): ", packageS[which( packageS %in% rownames(installed.packages()) == "FALSE")], " are not installed!")
  }
  
  
  library(csaw)
  library(edgeR)
  library(GenomicAlignments)
  library(data.table)
  
  #################################################################################################################
  
  ## Enter specified working directory:
  message("")
  
  if (NAME == "") stop("Please enter a sample name!")
  
  if (WD != "~/"){
    
    if (dir.exists(WD)){
      message("Enter working directory ", WD)
      setwd(WD)
    }
    
    if (!dir.exists(WD)){
      message("Creating working directory ", WD)
      dir.create(WD, showWarnings = T)
      setwd(WD)
    }
    
  }
  
  #################################################################################################################
  ## message if no name given in function
  if ( (NAME == "") || (CONTRASTS == "") || (GROUP == "") ){
    stop("Please enter Name/Contrasts/Group for the analysis!")
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
  #  
  #  out <- getPESizes(BAMS[1])out
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
  
  assign( paste(NAME, "_regionCounts", sep=""),
          data, envir = .GlobalEnv)
  
  #################################################################################################################
  ## Normalize 
  message("Calculating normalization factors")
  
  ## Save counts and put in colnames:
  strReverse <- function(x){
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
  }
  file.names <- as.character(strReverse(sapply(strReverse(BAMS), function(x) strsplit(x, split = "\\/")[[1]][1])))
  
  RAWcounts <- data.frame(assay(data))
  colnames(RAWcounts) <- file.names
  
  ## save normalized counts based on peak counts
  CPMcountsPeaks <- data.frame( calculateCPM(object = normFactors(data, se.out=data), use.norm.factors = T, log = F) )
  colnames(CPMcountsPeaks) <- file.names
  
  
  ## save normalized counts based on large bins
  binned <- windowCounts(BAMS, bin=TRUE, width=10000, param=PARAM)
  
  assign( paste(NAME, "_binned", sep=""),
          binned, envir = .GlobalEnv)
  
  CPMcountsBins  <- data.frame( calculateCPM(object = normFactors(binned, se.out=data), use.norm.factors = T, log = F) )
  colnames(CPMcountsBins) <- file.names
  
  ## depending on the settings of NORM use either peaks or bins for downstream analysis
  if (NORM == "peaks") {
    data <- normFactors(data, se.out=data)
  }
  
  if (NORM == "largebins") {
    data <- normFactors(binned, se.out=data)
  }
  
  ## Save counts:
  assign( paste(NAME, "_countsRAW", sep=""),
          RAWcounts, envir = .GlobalEnv)
  
  assign( paste(NAME, "_countsCPM_peaks", sep=""),
          CPMcountsPeaks, envir = .GlobalEnv)
  
  assign( paste(NAME, "_countsCPM_largebins", sep=""),
          CPMcountsBins, envir = .GlobalEnv)

  #################################################################################################################
  
  ## edgeR part:
  message("Running edgeR")
  
  ## set up DGElist:
  y <- asDGEList(data)
  y$samples$group <- GROUP
  
  ## design matrix:  
  design <- model.matrix(~ 0 + group, data=y$samples)
  colnames(design) <- gsub("group", "", colnames(design))
  assign( paste(NAME, "_design", sep=""),
          design, envir = .GlobalEnv)
  
  ## Model fitting:
  ## estimate dispersion:
  message("Estimating dispersion")
  y <- estimateDisp(y, design)
  assign( paste(NAME, "_estimateDisp", sep=""),
          y, envir = .GlobalEnv)
  
  message("Fitting GLM")
  fit <- glmQLFit(y, design, robust=TRUE)
  assign( paste(NAME, "_fit", sep=""),
          fit, envir = .GlobalEnv)
    
  #################################################################################################################
  ## test all contrasts:
  message("Testing all contrasts")
  for (i in CONTRASTS){
    
   ## current contrast:
   current.results <- glmQLFTest(fit, contrast = makeContrasts(contrasts = i, levels=design))
   assign( paste(NAME, "_resultsRAW_", gsub("-", "_",i), sep=""),
           current.results, envir = .GlobalEnv)
      
   ## fdr:
   current.out <- topTags(current.results, n=Inf, adjust.method="BH", sort.by="none")
   current.out <- data.frame( data.frame(rowRanges(data))[,1:3], current.out$table)
   assign( paste(NAME, "_resultsFDR_", gsub("-", "_",i), sep=""),
           current.out, envir = .GlobalEnv)

   }
  
  #################################################################################################################
  ## Produce MA plots, using the average per replicate group:
  plotMA_custom <- function(COUNTS, MAIN = "", REPLACE.ZERO = 1){
    
    R=COUNTS[,1]
    G=COUNTS[,2]  
    
    R[which(R == 0)] <- REPLACE.ZERO
    G[which(G == 0)] <- REPLACE.ZERO
    
    ## get mean of counts and log2FC
    M <- log2(R/G)
    A <- 0.5*log2(R*G)
    
    ## Decide y-axis limits based on rounded quantiles
    YLIM <- c(floor(quantile(M, 0.0001, na.rm=T)), ceiling(quantile(M, 0.9999, na.rm=T)))
    
    ## using smoothScatter
    par(bty="n")
    smoothScatter(A, M, main = MAIN, bty="n",
                  xlab="mean of normalized counts",
                  ylab="log2FC", ylim=YLIM, nbin = 200)
    abline(h=0)  
    
    
    ## plot points beyond y-axis limits as triangles to the y-axis limits
    points(x = A[which(M < YLIM[1])], y = rep(YLIM[1], length(which(M < YLIM[1]))), pch=17, cex = 0.7, col="darkblue")
    points(x = A[which(M > YLIM[2])], y = rep(YLIM[2], length(which(M > YLIM[2]))), pch=17, cex = 0.7, col="darkblue")
    
  }
  
  Do_MAplot <- function(CPMs, SUFFIX){
    if (MAall == F){
      message("Printing MA plots averaged over replicates")
      
      
      if ("FALSE" %in% grepl("rep", colnames(CPMs))) message("Warning: Not all sample names end with _rep* -- skipping MAplot part!")
      
      ## If naming convention *_rep* is ok, proceed:
      if( length( grep("FALSE", grepl("rep", colnames(CPMs))) ) == 0){
        
        ## Average all replicates per condition:
        names.unique <- unique(sapply(strsplit(colnames(CPMs), split="_rep"), function(x) x[1]))
        
        ## Make all pairwise comparisons:
        tmp.combn <- combn(x = names.unique, m = 2)
        
        pdf(file = paste(NAME,"_MAplots_",SUFFIX,".pdf", sep="") , onefile=T, paper='A4') 
        par(mfrow=c(1,1), bty="n")
        for (i in 1:ncol(tmp.combn)){
          
          one <- CPMs[,grep(tmp.combn[1,i], colnames(CPMs))]
          two <- CPMs[,grep(tmp.combn[2,i], colnames(CPMs))]
          
          if (ncol(one) == 1) average.one <- one
          if (ncol(one) > 1) average.one <- rowMeans(one)
          if (ncol(two) == 1) average.two <- two
          if (ncol(two) > 1) average.two <- rowMeans(two)
          
          ## MA:
          tmp.df <- data.frame(average.one, average.two)
          colnames(tmp.df) <- tmp.combn[,i]
          
          
          plotMA_custom(COUNTS = tmp.df, MAIN = paste("Sample: ", paste(tmp.combn[,i], collapse=" vs "), sep=""))
          
        }; dev.off()
      }
    }
    
    if (MAall == T){
      
      message("Printing MA plots for all sample combinations")
      
      tmp.combn <- combn(x = colnames(CPMs), m = 2)
      
      pdf(file = paste(NAME,"_MAplots_",SUFFIX,".pdf", sep="") , onefile=T, paper='A4') 
      par(mfrow=c(1,1), bty="n")
      for (i in 1:ncol(tmp.combn)){
        
        one <- CPMs[,grep(tmp.combn[1,i], colnames(CPMs))]
        two <- CPMs[,grep(tmp.combn[2,i], colnames(CPMs))]
        
        ## MA:
        tmp.df <- data.frame(one, two)
        colnames(tmp.df) <- tmp.combn[,i]
        
        
        plotMA_custom(COUNTS = tmp.df, MAIN = paste("Sample: ", paste(tmp.combn[,i], collapse=" vs "), sep=""))
        
      }; dev.off()
    }
  }
  Do_MAplot(CPMs = CPMcountsBins, SUFFIX = "largebins")
  Do_MAplot(CPMs = CPMcountsPeaks, SUFFIX = "peakbased")
  
  
  #################################################################################################################
  message("[ENDED]: ", NAME)
  message("")
  
}
