## Function accepts raw count matrix as GRanges.
## Will then normalize using TMM, make MA-plots to check efficiency and do a basic PCA.
## If a GRanges with counts from genome-wide bins is given uses this for normalization factor calculation.

###########################################################################################################

csaw_counts2cpm <- 
  function(Basename,                       ## the prefix for everything written to disk
           Counts.gr,                      ## GRanges with raw counts
           Bins.gr = NULL,                 ## if a count matrix is provided here, use for largebins norm.
           Colnames.Suffix = "_dedup.bam", ## either suffix or ""
           MakeMAplots = "groupwise",      ## c("all", "groupwise", "none")
           MakePCAplot = TRUE,             ## self explicatory
           WorkingDir = "./"               ## the directory to save MA-plots
  )
  {
  
  GetDate <- function(){ gsub("^20", "", format(Sys.Date(), "%Y%m%d")) }
  
  #########################################################################################################
  
  ## Check if required packages are installed:
  packageS <- c("csaw")
  
  if (length(grep("FALSE", (packageS %in% rownames(installed.packages())))) > 0){
    stop("==> Package(s): ", 
         packageS[which( packageS %in% rownames(installed.packages()) == "FALSE")], 
         " are not installed!")
  }
  
  library(csaw)
  
  #########################################################################################################
  
  ## Check parameters:
  
  if (class(Counts.gr) != "GRanges")         stop("==> Counts.gr must be GRanges")
  if (ncol(elementMetadata(Counts.gr)) == 0) stop("==> Counts.gr appears to have no count data")
  
  if (!is.numeric(as.matrix(elementMetadata(counts.gr)))) {
    stop("==> Counts.gr seems to contain non-numeric values")
  }
  
  
  ## This one checks if any unwanted columns made it to the count matrix in Counts.gr:
  if (sum(!grepl("_rep", colnames(elementMetadata(Counts.gr)))) > 0) {
    warning("==> Not all colnames in the count matrix of Counts.gr have <_rep> in it. Better check that!")
  }
  
  if (!is.null(Bins.gr)){
    
    if (class(Bins.gr) != "GRanges")   stop("==> Bins.gr must be GRanges")
    
    if (ncol(elementMetadata(Bins.gr)) == 0) stop("==> Bins.gr appears to have no count data")
    
    if (sum(!grepl("_rep", colnames(elementMetadata(Bins.gr)))) > 0) {
      warning("==> Not all colnames in the count matrix of Bins.gr have <_rep> in it. Better check that!")
    }
    
    if (!is.numeric(as.matrix(elementMetadata(Bins.gr)))) {
      stop("==> Bins.gr seems to contain non-numeric values")
    }
    
  }
  
  ## Check if MakeMAplots is one of the three supported values:
  if (sum(grepl(pattern = paste("^", MakeMAplots, "$", sep=""), 
                x = c("all", "groupwise", "none"))
  ) == 0){ stop("==> MakeMAplots must be one of (all, groupwise, none)")  }
  
  if (!is.logical(MakePCAplot)){
    warning("==> MakePCAplot can only be TRUE or FALSE, set to FALSE now.")
    MakePCAplot <- FALSE
  }
  
  #########################################################################################################
  
  ## If not exists, create directories for plots and lists:
  if (!dir.exists(WorkingDir)){
    message("Creating WorkingDir: -- ", normalizePath(WorkingDir))
  }
  
  if (!dir.exists(paste(WorkingDir, "Lists",sep=""))){
    dir.create(paste(WorkingDir, "Lists",sep=""), showWarnings = T)
  }
  
  if (!dir.exists(paste(WorkingDir, "Plots",sep=""))){
    dir.create(paste(WorkingDir, "Plots",sep=""), showWarnings = T)
  }
  
  #########################################################################################################
  
  message("[Started]:","\t", Basename)
  
  ## Remove suffix from colnames:
  colnames(elementMetadata(Counts.gr)) <- gsub(Colnames.Suffix, 
                                               "", 
                                               colnames(elementMetadata(Counts.gr)))
  
  ## Create SummarizedExperiments object from input count matrix (which is GRanges):
  data.se <- SummarizedExperiment(
    list(counts=as.matrix(elementMetadata(Counts.gr))))
  
  ## add total depth and ranges:
  data.se$totals     <- colSums(as.matrix(elementMetadata(Counts.gr)))
  rowRanges(data.se) <- Counts.gr[,0]
  
  #########################################################################################################
  
  ## If no counts for large bins are provided use TMM based on peak counts:
  if (is.null(Bins.gr) == TRUE){
    message("Normalizing data based on peak counts")
    data.se <- normFactors(data.se, se.out=data.se)
  }
  
  ## Else, use bins:
  if (is.null(Bins.gr) == FALSE){
    
    message("Normalizing data based on large bins")
    
    colnames(elementMetadata(Bins.gr)) <- gsub(Colnames.Suffix, 
                                               "", 
                                               colnames(elementMetadata(Bins.gr)))
    
    bins.se <- SummarizedExperiment(
      list(counts=as.matrix(elementMetadata(Bins.gr))))
    
    bins.se$totals     <- colSums(as.matrix(elementMetadata(Bins.gr)))
    
    nF <- normFactors(bins.se, se.out=FALSE)
    
    data.se$norm.factors <- nF
  
  }
  
  #########################################################################################################
  
  ## Calculate per-millions:
  message("Calculating CPMs")
  CPM.counts.gr <- Counts.gr[,0]
  
  elementMetadata(CPM.counts.gr) <- as.matrix(
    calculateCPM(object = data.se, use.norm.factors = T, log = F))
  
  ## Save cpm to disk:
  message("Writing CPMs to disk")
  write.table(x = data.frame(CPM.counts.gr)[,-c(4,5)], 
              file = paste(WorkingDir, "Lists/", 
                           GetDate(),"_", Basename, "_normCounts.tsv",
                           sep=""),
              quote = F, 
              col.names = T, 
              row.names = F,
              sep="\t")
  
  #########################################################################################################
  
  ## Main Function to produce MA-plots:
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
    points(x = A[which(M < YLIM[1])], y = rep(YLIM[1], 
                                              length(which(M < YLIM[1]))), pch=17, cex = 0.7, col="darkblue")
    points(x = A[which(M > YLIM[2])], y = rep(YLIM[2], 
                                              length(which(M > YLIM[2]))), pch=17, cex = 0.7, col="darkblue")
    
  }
  
  #########################################################################################################
  
  Do_MAplot <- function(CPMs){
    
    ## MA-plots averaged over groups like *_rep*{Colnames.Suffix}, everything before _rep is a group:
    if (MakeMAplots == "groupwise"){
      message("Printing MA-plots averaged over replicate groups")
      
      if ("FALSE" %in% grepl("rep", colnames(CPMs))) {
        warning("Not all sample names end with _rep* -- skipping MAplot part!")
      }
      
      ## If naming convention *_rep* is ok, proceed:
      if( length( grep("FALSE", grepl("rep", colnames(CPMs))) ) == 0){
        
        ## Average all replicates per condition:
        names.unique <- unique(sapply(strsplit(colnames(CPMs), split="_rep"), function(x) x[1]))
        
        ## Make all pairwise comparisons:
        tmp.combn <- combn(x = names.unique, m = 2)
        
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
      }
    }
    
    if (MakeMAplots == "all"){
      
      message("Printing MA-plots for all sample combinations")
      
      tmp.combn <- combn(x = colnames(CPMs), m = 2)
      
      for (i in 1:ncol(tmp.combn)){
        
        one <- CPMs[,grep(tmp.combn[1,i], colnames(CPMs))]
        two <- CPMs[,grep(tmp.combn[2,i], colnames(CPMs))]
        
        ## MA:
        tmp.df <- data.frame(one, two)
        colnames(tmp.df) <- tmp.combn[,i]
        
      }
    }
  }
  
  ## Call MA-plot function and append a suffix to the final pdf that indicates which norm method was used:
  if (!is.null(Bins.gr)) SUFFIX <- "largebins"  
  if (is.null(Bins.gr))  SUFFIX <- "peakbased"
  
  if (MakeMAplots != "none"){
    
    par(mfrow=c(1,1), bty="n")
    pdf(paper = "a4", file = paste(WorkingDir,"./Plots/", GetDate(), "_MAplots_", SUFFIX, ".pdf", sep=""))
    
    Do_MAplot(CPMs = data.frame(elementMetadata(CPM.counts.gr)))
    
    dev.off()
  
  }
  
  #########################################################################################################
  ## Plot a PCA (for this transform cpm to rangedsummarizedexperiment and add coldata):
  if (MakePCAplot == TRUE){
    
    message("Plotting PCA")
    
    rse <- SummarizedExperiment(
      list(counts=as.matrix(elementMetadata(CPM.counts.gr))))
    
    rowRanges(rse) <- CPM.counts.gr[,0]
    rse <- DESeqTransform(rse)
    rownames(colData(rse)) <- colnames(rse)
    colData(rse)$group <- as.factor(sapply(strsplit(colnames(rse), split="_rep"), function(x)x[1]))
    
    pdf(paper = "a4", file = paste(WorkingDir, "/Plots/", GetDate(), "_PCA_top2k.pdf", sep=""))
    print(plotPCA(object = rse, intgroup="group", 2000))
    dev.off()
  }
    
  #########################################################################################################
  
  message("[Ended]:","\t", Basename)
  
  #########################################################################################################
  
}
