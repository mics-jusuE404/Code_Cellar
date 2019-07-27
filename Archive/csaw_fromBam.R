## Some code to get count matrix from bam files with csaw:

################################################################################################################################################

library(csaw)
library(GenomicAlignments)
library(data.table)

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
  
  ## Also save to disk:
  if (!dir.exists("./Lists")) dir.create("./Lists")
  
  GetDate <- function(){ gsub("^20", "", format(Sys.Date(), "%Y%m%d")) }
  
  write.table(data.frame(raw.gr), quote = F, col.names = T, row.names = F, sep="\t", 
              file = paste("./Lists/", GetDate(), "_", NAME, "_countsRAW.tsv", sep=""))
  write.table(data.frame(cpm.gr), quote = F, col.names = T, row.names = F, sep="\t", 
              file = paste("./Lists/", GetDate(), "_", NAME, "_countsCPM.tsv", sep=""))
