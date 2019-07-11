## Template for differential count analysis (ChIP-seq, ATAC-seq) with csaw/edgeR.
## Three functions are used:
## 1) run_csaw_peakBased to count reads over peaks, normalization and producing MA-plots,
## 2) edgeR_FitGLM to estimate dispersion and fit the glm based on the experimental design,
## 3) edgeR_TestContrasts to test the specified contrasts, optionally against a certain FC with glmTreat()

################################################################################################################################################
################################################################################################################################################

library(csaw)
library(edgeR)
library(GenomicAlignments)
library(data.table)
library(statmod)

################################################################################################################################################
################################################################################################################################################
## Function counts/normalizes reads and fits the GLM. In a separate function one then tests the contrasts:

run_csaw_peakBased <- function(NAME,                         ## the name assigned to this analysis
                               SUMMITS,                      ## path to a BED-like file with summits (or intervals)
                               BLACKLIST,                    ## path to a blacklist excluded for analysis
                               WIDTH = 200,                  ## window size to extend summits
                               BAMS,                         ## full path to BAM files
                               FRAGLEN = "",                 ## length to extend reads to fragments (use 1 for ATAC-seq)
                               PAIRED = F,                   ## if using paired-end mode (so far not implemented)
                               NORM = "peaks",               ## "peaks" or "largebins" for TMM-normalization
                               FILTER_aveLogCPM = c(-1) ,    ## apply aveLogCPM filter default -1, set to "" to deactivate
                               CORES= c(detectCores()-1),    ## number of cores for read counting, default is all but one
                               plotMAall = F,                ## plot no (none), all possible (T) or group-wise (N) MA-plots
                               PLOTDIR,                      ## the directory to save MA-plots
                               FILENAMES = ""                ## file names to be used in MA plots (can be helpful if BAM files contain special chars)
                               ){              ## whether to calculate and save normalized counts based on large bins
  
  ## Check if required packages are installed:
  packageS <- c("csaw", "statmod", "edgeR", "GenomicAlignments", "data.table")
  if (length(grep("FALSE", (packageS %in% rownames(installed.packages())))) > 0){
    stop("Package(s): ", packageS[which( packageS %in% rownames(installed.packages()) == "FALSE")], " are not installed!")
  }
  
  if (PAIRED == T) stop("Paired-end mode not yet implemented")
  
  if (FILTER_aveLogCPM != "" && is.numeric(FILTER_aveLogCPM) == F) stop("FILTER_aveLogCPM is not numeric")
  
  ##############################################################################################################################
  ## Enter specified working directory (or create if not existing):
  message("")
  message("[STARTED]: ", NAME)
  
  if (!dir.exists(PLOTDIR)){
    message("Creating plotting directory ", PLOTDIR)
    dir.create(PLOTDIR, showWarnings = T)
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
  
  
  if (FILTER_aveLogCPM != ""){
    
    message("Filtering regions for aveLogCPM of ", FILTER_aveLogCPM)
    keep <- aveLogCPM(asDGEList(data)) >= FILTER_aveLogCPM
    data <- data[keep,]
    
  }

  #################################################################################################################
  ## Function returns the file name (essentially everything in full path after the last "/"):
  strReverse <- function(x){
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
  }
  
  if (FILENAMES[1] == ""){
    file.names <- as.character(strReverse(sapply(strReverse(BAMS), function(x) strsplit(x, split = "\\/")[[1]][1])))
  }
  if (FILENAMES[1] != ""){
    file.names <- FILENAMES
  }
  #################################################################################################################
  
  ## Normalize according to NORM parameter (peaks or bins as reference, default is peaks)
  if (NORM == "peaks") {
    message("Normalize counts based on peak counts")
    data <- normFactors(data, se.out=data)
  }
  
  if (NORM == "largebins") {
    message("Normalize counts based on bin counts")
    binned <- windowCounts(BAMS, bin=TRUE, width=10000, BPPARAM = MulticoreParam(workers = CORES), param = PARAM)
    data <- normFactors(binned, se.out=data)
  }
  
  ## store raw counts:
  RAWcounts <- data.frame(assay(data))
  colnames(RAWcounts) <- gsub("_dedup.bam", "", file.names)
  
  ## store CPMs based on TMM with peak counts:
  CPMcounts <- data.frame( calculateCPM(object = normFactors(data, se.out=data), use.norm.factors = T, log = F) )
  colnames(CPMcounts) <- gsub("_dedup.bam", "", file.names)
  
  ## Save raw and CPM counts as GRanges:
  raw.gr <- makeGRangesFromDataFrame(df = cbind( rowRanges(data), RAWcounts),
                                     keep.extra.columns = T)
  
  cpm.gr <- makeGRangesFromDataFrame(df = cbind( rowRanges(data), CPMcounts),
                                     keep.extra.columns = T)

  assign( paste(NAME, "_countsRAW.gr", sep=""),
          raw.gr, envir = .GlobalEnv)
  
  assign( paste(NAME, "_countsCPM.gr", sep=""),
          cpm.gr, envir = .GlobalEnv)
  
  assign( paste(NAME, "_data", sep=""),
          data, envir = .GlobalEnv)
  
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
        
        pdf(file = paste(PLOTDIR, paste(NAME,"_MAplots_",SUFFIX,".pdf", sep=""), sep="/") , onefile=T, paper='A4') 
        par(mfrow=c(1,1), bty="n")
        for (i in 1:ncol(tmp.combn)){
          
          one <- CPMs[,grep(tmp.combn[1,i], colnames(CPMs))]
          two <- CPMs[,grep(tmp.combn[2,i], colnames(CPMs))]
          
          if (ncol(data.frame(one)) == 1) average.one <- one
          if (ncol(data.frame(one)) > 1) average.one <- rowMeans(one)
          if (ncol(data.frame(two)) == 1) average.two <- two
          if (ncol(data.frame(two)) > 1) average.two <- rowMeans(two)
          
          ## MA:
          tmp.df <- data.frame(average.one, average.two)
          colnames(tmp.df) <- tmp.combn[,i]
          
          
          plotMA_custom(COUNTS = tmp.df, MAIN = paste("Sample: ", paste(tmp.combn[,i], collapse=" vs "), sep=""))
          
        }
        dev.off()
      }
    }
    
    if (plotMAall == T){
      
      message("Printing MA plots for all sample combinations")
      
      tmp.combn <- combn(x = colnames(CPMs), m = 2)
      
      pdf(file = paste(PLOTDIR, paste(NAME,"_MAplots_",SUFFIX,".pdf", sep=""), sep="/") , onefile=T, paper='A4') 
      par(mfrow=c(1,1), bty="n")
      for (i in 1:ncol(tmp.combn)){
        
        one <- CPMs[,grep(tmp.combn[1,i], colnames(CPMs))]
        two <- CPMs[,grep(tmp.combn[2,i], colnames(CPMs))]
        
        ## MA:
        tmp.df <- data.frame(one, two)
        colnames(tmp.df) <- tmp.combn[,i]
        
        
        plotMA_custom(COUNTS = tmp.df, MAIN = paste("Sample: ", paste(tmp.combn[,i], collapse=" vs "), sep=""))
        
      }
      dev.off()
    }
  }
  
  if (plotMAall != "none") {
    
    if (NORM == "peakbased") suf <- "peaks"
    if (NORM == "largebins") suf <- "largebins"
    
    Do_MAplot(CPMs = CPMcounts, SUFFIX = suf)
    
  }
  
  
  #################################################################################################################
  message("[ENDED]: ", NAME)
  message("")
  
}

################################################################################################################################################
################################################################################################################################################
## Function to estimate dispersion and fitting GLM:
edgeR_FitGLM <- function(DATA,      ## the *_regionCounts object from the above function
                         DESIGN,    ## the experimental design
                         NAME){     ## name as above
  
  ## set up DGElist:
  y <- asDGEList(DATA)
  
  ## Dispersion estimates:
  message("Estimating dispersion")
  y <- estimateDisp(y, DESIGN)
  assign( paste(NAME, "_estimateDisp", sep=""),
          y, envir = .GlobalEnv)
  
  message("Fitting GLM")
  fit <- glmQLFit(y, DESIGN, robust=TRUE)
  assign( paste(NAME, "_glmQLFit", sep=""),
          fit, envir = .GlobalEnv)
  
}

################################################################################################################################################
################################################################################################################################################

## Function to test the contrasts:
edgeR_TestContrasts <- function(CONTRASTS,       ## the output of makeContrasts()
                                FIT,             ## the output of glmQLFit() from the above function
                                GLMTREAT.FC="",  ## if numeric input test against that FC (FC not logFC) using glmTreat
                                RANGES,          ## the *_regionCounts object from above to get the genomic ranges
                                NAME){           ## name for assign()
  
  message("Testing all contrasts (total of ", dim(CONTRASTS)[2],")")
  
  for (i in seq(1,dim(CONTRASTS)[2])){
    
    ## current contrast:
    if (GLMTREAT.FC == ""){
      message("Null hypothesis is FC=0")
      current.results <- glmQLFTest(FIT, contrast = CONTRASTS[,i])
    }
    if (GLMTREAT.FC != "" && is.numeric(GLMTREAT.FC)){
      message("Null hypothesis is FC = ",GLMTREAT.FC, " for ", attr(CONTRASTS, "dimnames")$Contrasts[i])
      current.results <- glmTreat(glmfit = FIT, contrast = CONTRASTS[,i], lfc = log2(GLMTREAT.FC))
    }
    
    ## Save the FDR-adjusted TT:
    current.out <- topTags(current.results, n=Inf, adjust.method="BH", sort.by="none")
    current.out <- data.frame( data.frame(rowRanges(RANGES))[,1:3], current.out$table)
    assign( paste(NAME, "_topTags_", gsub("-", "_", attr(CONTRASTS, "dimnames")$Contrasts[i]), sep=""),
            current.out, envir = .GlobalEnv)
    rm(current.results)
    
  }
}

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

#### EXAMPLES:

## List BAM files:
#tmp.bam  <- list.files("~/path/to/bam/directory/", full.names = T, pattern = "\\.bam$")

## Extract the names only:
#tmp.name <- gsub(".bam", "", sapply(strsplit(tmp.bam, split = "bam/"), function(x) x[2]))#

## Define coldata indicating groups, genotypes, whatever defines the samples:
#tmp.coldata <- data.frame(NAME      = tmp.name,
#                           GENOTYPE  = sapply(strsplit(tmp.name, split="_"), function(x) x[1]),
#                           CELLTYPE  = sapply(strsplit(tmp.name, split="_"), function(x) x[2]),
#                           FACTORIAL = sort(rep(c("A", "B", "C", "D"), 2)))

## Design matrix for edgeR without intercept to have full control over contrasts:
#tmp.design <- model.matrix(~ 0 + FACTORIAL, data=tmp.coldata)#

## Specify contrasts for all comparisons:
#tmp.contrasts <- makeContrasts(Contr.Interaction = (FACTORIALA-FACTORIALB) - (FACTORIALC-FACTORIALD),
                                #Contr.Average     = (FACTORIALA+FACTORIALC)/2 - (FACTORIALB+FACTORIALD)/2,
                                #levels = tmp.design)

## Count and normalize reads over specified peaks:
#run_csaw_peakBased(NAME = "test", SUMMITS = "~/IMTB/Our_Data/Fischer2019/ATAC-seq/peaks/ATACseq_combinedCall_summits.bed", 
#                   BLACKLIST = "/Volumes/Rumpelkammer/Genomes/mm10/mm10_consensusBL.bed", WIDTH = 200, 
#                   BAMS = tmp.bam, FRAGLEN = 1, PAIRED = F, NORM = "peaks",
#                   CORES = 16, plotMAall = F, PLOTDIR = "~/testdir")

## Estimate disp and fit GLM:
#edgeR_FitGLM(DATA = test_regionCounts, DESIGN = tmp.design, NAME = "test")

## Test contrasts optionally against a minimum FC:
#edgeR_TestContrasts(CONTRASTS = tmp.contrasts, FIT = test_glmQLFit, GLMTREAT.FC = 2, RANGES = test_regionCounts, NAME = "test")
