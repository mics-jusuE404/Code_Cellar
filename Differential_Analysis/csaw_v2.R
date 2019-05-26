## Template for differential count analysis of DNA-seq data such as ATAC-seq or ChIP-seq.
## In its current state, it uses csaw (with user-provided peak summit positions, extended to user-defined size)
## to count reads and perform edgeR-based differential analysis.
## It will output all relevant intermediate results and CPMs plus MA-plots either comparing
## every possible combination of samples after normalization or averaged based on defined groups.

run_csaw_peakBased <- function(NAME,                         ## the name assigned to this analysis
                               SUMMITS,                      ## path to a BED-like file with summits (or intervals)
                               BLACKLIST,                    ## path to a blacklist excluded for analysis
                               WIDTH = 200,                  ## window size to extend summits
                               BAMS,                         ## full path to BAM files
                               FRAGLEN = "",                 ## length to extend reads to fragments (use 1 for ATAC-seq)
                               PAIRED = F,                   ## if using paired-end mode (so far not implemented)
                               NORM = "peaks",               ## "peaks" or "largebins" for TMM-normalization
                               DESIGN,                       ## design fulfillign edgeR requirements
                               CONTRASTS,                    ## contrasts to test
                               CORES= c(detectCores()-1),    ## number of cores for read counting, default is all but one
                               plotMAall = F,                ## plot no (none), all possible (T) or group-wise (N) MA-plots
                               WORKINGDIR,                   ## the directory to save MA-plots
                               saveBinned = F){              ## whether to calculate and save normalized counts based on large bins
  
  ## Check if required packages are installed:
  packageS <- c("csaw", "statmod", "edgeR", "GenomicAlignments", "data.table")
  if (length(grep("FALSE", (packageS %in% rownames(installed.packages())))) > 0){
    stop("Package(s): ", packageS[which( packageS %in% rownames(installed.packages()) == "FALSE")], " are not installed!")
  }
  
  library(csaw)
  library(edgeR)
  library(GenomicAlignments)
  library(data.table)
  library(statmod)
  
  if (PAIRED == T) stop("Paired-end mode not yet implemented")
  
  if (NORM == "largebins") saveBinned = T
  ##############################################################################################################################
  ## Enter specified working directory (or create if not existing):
  message("")
  message("[STARTED]: ", NAME)
  
  if (dir.exists(WORKINGDIR)){
    message("Enter working directory ", WORKINGDIR)
    setwd(WORKINGDIR)
  }
  
  if (!dir.exists(WORKINGDIR)){
    message("Creating working directory ", WORKINGDIR)
    dir.create(WORKINGDIR, showWarnings = T)
    setWORKINGDIR(WORKINGDIR)
  }

  #################################################################################################################
  ## Extend peak summits to granges:
  tmp.summits <- fread(SUMMITS, data.table = F, header = F)
  tmp.gr <- GRanges(seqnames = tmp.summits$V1,
                    ranges = IRanges(start = tmp.summits$V2+1, end = tmp.summits$V3))
  
  ## Resize and merge in case of overlaps:
  message("Extending peak summits to ", WIDTH, "bp")
  peaks_resized.gr <- GenomicRanges::reduce(GenomicRanges::resize(tmp.gr, width = WIDTH, fix = "center"))
  
  #################################################################################################################
  ## Fragment length:
  #if (PAIRED == T) {
  #  
  #  out <- getPESizes(BAMS[1])out
  #}
  
  if (PAIRED == F) PSE = "none"
  
  ## PARAM, discarding blacklisted regions
  PARAM=readParam(discard = makeGRangesFromDataFrame(df = fread(BLACKLIST, data.table=F), 
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
  message("Counting reads per peaks for ", length(BAMS), " BAM files using ", CORES, " workers")
  data <- regionCounts(BAMS, 
                       ext = FRAGLEN, 
                       param = PARAM,
                       regions = peaks_resized.gr,
                       BPPARAM = MulticoreParam(workers = CORES))
  
  assign( paste(NAME, "_regionCounts", sep=""),
          data, envir = .GlobalEnv)
  
  #################################################################################################################
  ## Function returns the file name (essentially everything in full path after the last "/"):
  strReverse <- function(x){
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
  }
  file.names <- as.character(strReverse(sapply(strReverse(BAMS), function(x) strsplit(x, split = "\\/")[[1]][1])))
  
  if (saveBinned == T){
  
    ## store CPMs based on TMM with largebin counts:
    message("Counting reads per bins for ", length(BAMS), " BAM files using ", CORES, " workers")
    binned <- windowCounts(BAMS, bin=TRUE, width=10000, param=PARAM)
 
    CPMcountsBins  <- data.frame( calculateCPM(object = normFactors(binned, se.out=data), use.norm.factors = T, log = F) )
    colnames(CPMcountsBins) <- file.names
  }
  
  message("Calculating normalization factors and storing counts")
  ## store raw counts:
  RAWcounts <- data.frame(assay(data))
  colnames(RAWcounts) <- file.names
  
  ## store CPMs based on TMM with peak counts:
  CPMcountsPeaks <- data.frame( calculateCPM(object = normFactors(data, se.out=data), use.norm.factors = T, log = F) )
  colnames(CPMcountsPeaks) <- file.names
  
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
  
  if (saveBinned == T){
    assign( paste(NAME, "_countsCPM_largebins", sep=""),
            CPMcountsBins, envir = .GlobalEnv)
  }
  
  #################################################################################################################
  ## Differential analysis with edgeR:
  message("Differential analysis with edgeR")
  
  ## set up DGElist:
  y <- asDGEList(data)
  
  ## Dispersion estimates:
  message("Estimating dispersion")
  y <- estimateDisp(y, DESIGN)
  assign( paste(NAME, "_estimateDisp", sep=""),
          y, envir = .GlobalEnv)
  
  message("Fitting GLM")
  fit <- glmQLFit(y, DESIGN, robust=TRUE)
  assign( paste(NAME, "_glmQLFit", sep=""),
          fit, envir = .GlobalEnv)
  
  #################################################################################################################
  ## Test all specified contrasts:
  message("Testing all contrasts (total of ", dim(CONTRASTS)[2],")")
  
  for (i in seq(1,dim(CONTRASTS)[2])){
    
    ## current contrast:
    current.results <- glmQLFTest(fit, contrast = CONTRASTS[,i])
    
    ## Save the raw result:
    assign( paste(NAME, "_glmQLFTest_", gsub("-", "_", attr(CONTRASTS, "dimnames")$Contrasts[i]), sep=""),
            current.results, envir = .GlobalEnv)
    
    ## Save the FDR-adjusted TT:
    current.out <- topTags(current.results, n=Inf, adjust.method="BH", sort.by="none")
    current.out <- data.frame( data.frame(rowRanges(data))[,1:3], current.out$table)
    assign( paste(NAME, "_topTags_", gsub("-", "_", attr(CONTRASTS, "dimnames")$Contrasts[i]), sep=""),
            current.out, envir = .GlobalEnv)
    rm(current.results)
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
    
    if (plotMAall == F){
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
    
    if (plotMAall == T){
      
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
  
  if (plotMAall != "none") {
    
    if (saveBinned == T) Do_MAplot(CPMs = CPMcountsBins, SUFFIX = "largebins")
    Do_MAplot(CPMs = CPMcountsPeaks, SUFFIX = "peakbased")
    
  }
  
  
  #################################################################################################################
  message("[ENDED]: ", NAME)
  message("")
  
}

## Example:
tmp.bam1  <- list.files("~/IMTB/Fischer2019/ATAC-seq/bam", full.names = T, pattern = "\\.bam$")

tmp.name1 <- gsub(".bam", "", sapply(strsplit(tmp.bam1, split = "bam/"), function(x) x[2]))

tmp.coldata1 <- data.frame(NAME      = tmp.name1,
                           GENOTYPE  = sapply(strsplit(tmp.name1, split="_"), function(x) x[1]),
                           CELLTYPE  = sapply(strsplit(tmp.name1, split="_"), function(x) x[2]),
                           FACTORIAL = sort(rep(c("A", "B", "C", "D"), 2)))

tmp.design1 <- model.matrix(~ 0 + FACTORIAL, data=tmp.coldata1)

tmp.contrasts1 <- makeContrasts(Contr.Interaction = (FACTORIALA-FACTORIALB) - (FACTORIALC-FACTORIALD),
                                Contr.Average     = (FACTORIALA+FACTORIALC)/2 - (FACTORIALB+FACTORIALD)/2,
                                levels = tmp.design1)

run_csaw_peakBased(NAME = "test", SUMMITS = "~/IMTB/Fischer2019/ATAC-seq/peaks/ATACseq_combinedCall_summits.bed", 
                   BLACKLIST = "/Volumes/Rumpelkammer/Genomes/mm10/mm10_consensusBL.bed", WIDTH = 200, 
                   BAMS = tmp.bam1, FRAGLEN = 1, PAIRED = F, NORM = "peaks", DESIGN = tmp.design1,
                   CONTRASTS = tmp.contrasts1, CORES = 16, plotMAall = F, WORKINGDIR = "~/testdir")
