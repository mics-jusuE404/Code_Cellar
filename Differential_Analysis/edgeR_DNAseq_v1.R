## Differential analysis for DNA-seq experiments (ChIP-seq, ATAC-seq) starting from a count matrix.
## Performs normalization, plots exploratory MA-plots and PCA, performs diff. analysis, saves everything
## to disk

###########################################################################################################

csaw_counts2TopTags <- 
  function(Basename,                       ## the prefix for everything written to disk
           Counts.gr,                      ## GRanges with raw counts
           Bins.gr = NULL,                 ## if a count matrix is provided here, use for largebins norm.
           Colnames.Suffix = "_dedup.bam", ## either suffix or ""
           MakeMAplots = "groupwise",      ## c("all", "groupwise", "none")
           MakePCAplot = TRUE,             ## PCA plot based of DESeq2::plotPCA
           WriteCpmToDisk = TRUE,          ##
           WorkingDir = "./",              ## the directory to save MA-plots
           Design = NULL,                  ## experimental design for diff. analysis
           FiltByExpr = FALSE,             ## whether or not to use that filter
           Contrasts = NULL,               ## contrasts to test
           Thresh.glmTreat = NULL,         ## minFC to test against
           Save.Image = TRUE,              ## whether to save Rdata to ./WorkingDir/R
           ImageName = "tmp"               ## name for the image to be saved
  )
  {
    
    GetDate <- function(){ gsub("^20", "", format(Sys.Date(), "%Y%m%d")) }
    
    #########################################################################################################
    
    ## Check if required packages are installed:
    packageS <- c("csaw", "edgeR", "statmod", "DESeq2")
    
    if (length(grep("FALSE", (packageS %in% rownames(installed.packages())))) > 0){
      stop(call. = FALSE, "==> Package(s): ", 
           packageS[which( packageS %in% rownames(installed.packages()) == "FALSE")], 
           " are not installed!")
    }
    
    library(csaw)
    library(edgeR)
    library(statmod)
    library(DESeq2)
    
    #########################################################################################################
    
    ## Check parameters:
    
    if (class(Counts.gr) != "GRanges")         stop(call. = FALSE, "==> Counts.gr must be GRanges")
    if (ncol(elementMetadata(Counts.gr)) == 0) stop(call. = FALSE, "==> Counts.gr appears to have no count data")
    
    if (!is.numeric(as.matrix(elementMetadata(Counts.gr)))) {
      stop(call. = FALSE, "==> Counts.gr seems to contain non-numeric values")
    }
    
    
    ## This one checks if any unwanted columns made it to the count matrix in Counts.gr:
    if (sum(!grepl("_rep", colnames(elementMetadata(Counts.gr)))) > 0) {
      warning("==> Not all colnames in the count matrix of Counts.gr have <_rep> in it. Better check that!")
    }
    
    if (!is.null(Bins.gr)){
      
      if (class(Bins.gr) != "GRanges")   stop(call. = FALSE, "==> Bins.gr must be GRanges")
      
      if (ncol(elementMetadata(Bins.gr)) == 0) stop(call. = FALSE, "==> Bins.gr appears to have no count data")
      
      if (sum(!grepl("_rep", colnames(elementMetadata(Bins.gr)))) > 0) {
        warning("==> Not all colnames in the count matrix of Bins.gr have <_rep> in it. Better check that!")
      }
      
      if (!is.numeric(as.matrix(elementMetadata(Bins.gr)))) {
        stop(call. = FALSE, "==> Bins.gr seems to contain non-numeric values")
      }
      
    }
    
    ## Check if MakeMAplots is one of the three supported values:
    if (sum(grepl(pattern = paste("^", MakeMAplots, "$", sep=""), 
                  x = c("all", "groupwise", "none"))
    ) == 0){ stop(call. = FALSE, "==> MakeMAplots must be one of (all, groupwise, none)")  }
    
    if (!is.logical(MakePCAplot)){
      warning("==> MakePCAplot can only be TRUE or FALSE, set to FALSE now.")
      MakePCAplot <- FALSE
    }
    
    if (!is.null(Thresh.glmTreat) & !is.numeric(Thresh.glmTreat)){
      stop("==> Thresh.glmTreat is not numeric")
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
    
    if (!dir.exists(paste(WorkingDir, "R",sep=""))){
      dir.create(paste(WorkingDir, "R",sep=""), showWarnings = T)
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
    CPM.Counts.gr <- Counts.gr[,0]
    
    elementMetadata(CPM.Counts.gr) <- as.matrix(
      calculateCPM(object = data.se, use.norm.factors = T, log = F))
    
    ## Save cpm to disk:
    if (WriteCpmToDisk == TRUE){
      message("Writing CPMs to disk")
      write.table(x = data.frame(CPM.Counts.gr)[,-c(4,5)], 
                  file = paste(WorkingDir, "Lists/", 
                               GetDate(),"_", Basename, "_normCounts.tsv",
                               sep=""),
                  quote = F, 
                  col.names = T, 
                  row.names = F,
                  sep="\t")
    }
    
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
            
            plotMA_custom(COUNTS = tmp.df, MAIN = paste("Sample: ", paste(tmp.combn[,i], collapse=" over "), sep=""))
            
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
      
      Do_MAplot(CPMs = data.frame(elementMetadata(CPM.Counts.gr)))
      
      dev.off()
      
    }
    
    #########################################################################################################
    ## Plot a PCA (for this transform cpm to rangedsummarizedexperiment and add coldata):
    if (MakePCAplot == TRUE){
      
      message("Plotting PCA")
      
      rse <- SummarizedExperiment(
        list(counts=as.matrix(elementMetadata(CPM.Counts.gr))))
      
      rowRanges(rse) <- CPM.Counts.gr[,0]
      rse <- DESeqTransform(rse)
      rownames(colData(rse)) <- colnames(rse)
      colData(rse)$group <- as.factor(sapply(strsplit(colnames(rse), split="_rep"), function(x)x[1]))
      
      pdf(paper = "a4", file = paste(WorkingDir, "/Plots/", GetDate(), "_PCA_top2k.pdf", sep=""))
      print(plotPCA(object = rse, intgroup="group", 2000))
      dev.off()
    }
    
    #########################################################################################################
    
    if (is.null(Design)) {
      message("No Design given therefore stop before differential analysis part")
      message("[Ended]:", "\t", Basename)
      #https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
      stop_quietly <- function() {
        opt <- options(show.error.messages = FALSE)
        on.exit(options(opt))
        stop()
      }; stop_quietly()
    }
    
    #########################################################################################################
    
    ## Start differential analysis:
    data.dge <- asDGEList(data.se)
    
    if (FiltByExpr == TRUE) {
      keep <- filterByExpr(data.dge, design = Design)
      data.se <- data.se[keep,]
      data.dge <- data.dge[keep, , keep.lib.sizes=FALSE]
      
    }
    
    ## dispersion estimates:
    message("Estimating dispersion")
    data.dge <- estimateDisp(data.dge, Design)
    
    ## glm data.fit:
    message("Fitting model")
    data.fit <- glmQLFit(data.dge, design = Design, robust = T)
    
    ## to envir:
    assign(paste(Basename, ".DGElist", sep=""),  data.dge, envir = .GlobalEnv)
    assign(paste(Basename, ".glmQLFit", sep=""), data.fit, envir = .GlobalEnv)
    
    #########################################################################################################
    
    ## before testing of contrasts starts define a function to plot MA-plots with signif regions colored:
    smoothScatter_TopTags <- function(TT, SIG.THRESHOLD=0.05, MAIN) {
      
      ## ylim based on quantiles:
      YLIM <- c(floor(quantile(TT$logFC, 0.0001, na.rm=T)), ceiling(quantile(TT$logFC, 0.9999, na.rm=T)))
      
      smoothScatter(x = TT$logCPM, y = TT$logFC,
                    xlab="average logCPM", 
                    ylab="log2FC", 
                    main=paste(MAIN,
                               " at ", 
                               "FDR = ", 
                               SIG.THRESHOLD*100, 
                               "%", sep=""), 
                    ylim=YLIM)
      
      ## count significant regions
      sig     <- TT[TT$FDR < SIG.THRESHOLD,]
      uppp    <- length(which(sig$logFC > 0))
      downnnn <- length(which(sig$logFC < 0))
      
      ## add significants regions in firebrick:
      points(x = sig$logCPM, y = sig$logFC, pch=20, col ="firebrick", cex=0.1)
      legend("topright", paste( paste("UP=",uppp, sep=""), 
                                paste("DOWN=", downnnn, sep=""), 
                                paste("minFC=", round(2^min(abs(sig$logFC)), digits=1), sep=""),
                                sep="\n"), bty="n") 
      
      
      COL<-"firebrick"
      
      ## add triangles for points beyond YLIM:
      lower.limit <- par("usr")[3]
      points(x = TT$logCPM[which(TT$logFC < lower.limit)], 
             y = rep( (lower.limit - lower.limit*0.025), length(which(TT$logFC < lower.limit))),
             pch=17, cex = 0.7, col=COL) 
      
      upper.limit <- par("usr")[4]
      points(x = TT$logCPM[which(TT$logFC > upper.limit)], 
             y = rep( (upper.limit - upper.limit*0.025), length(which(TT$logFC > upper.limit))),
             pch=17, cex = 0.7, col=COL) 
    }
    
    #########################################################################################################
    
    ## Test all Contrasts, optionally against a FC, and correct for multiple testing:
    for (i in seq(1,dim(Contrasts)[2])){
      
      ## Against 0:
      if (is.null(Thresh.glmTreat)){
        message("Null hypothesis is FC = 0")
        current.results <- glmQLFTest(data.fit, contrast = Contrasts[,i])
      }
      
      ## Against Thresh.glmTreat:
      if (!is.null(Thresh.glmTreat)){
        
        message("Null hypothesis is FC = ",
                Thresh.glmTreat, 
                " for ", 
                attr(Contrasts, "dimnames")$Contrasts[i])
        
        current.results <- glmTreat(glmfit = data.fit, 
                                    contrast = Contrasts[,i], 
                                    lfc = log2(Thresh.glmTreat))
        
      }
      
      ## Save the FDR-adjusted TT:
      current.out <- topTags(current.results, n=Inf, adjust.method="BH", sort.by="none")
      current.out <- data.frame( data.frame(rowRanges(data.se))[,1:3], current.out$table)
      
      if (is.null(Thresh.glmTreat))  thr <- 0
      if (!is.null(Thresh.glmTreat)) thr <- Thresh.glmTreat
      
      ## save as variable to envir indicating the null hypothesis in the name:
      tmp.name <- paste(Basename, 
                        "_topTags_", 
                        gsub("-", "_", 
                             attr(Contrasts, "dimnames")$Contrasts[i]),
                        "_null_",
                        thr,
                        sep="")
      
      ## as variable to envir:
      assign( tmp.name,
              current.out,
              envir = .GlobalEnv)
      
      ## to disk:
      write.table(current.out, sep="\t", col.names = T, row.names = F, quote = F,
                  file=paste(WorkingDir, "/Lists/", GetDate(), "_", tmp.name, ".tsv", sep=""))
      
      write.table(current.out[current.out$FDR < 0.05,], sep="\t", col.names = T, row.names = F, quote = F,
                  file=paste(WorkingDir, "/Lists/", GetDate(), "_", tmp.name, "_FDR5perc", ".tsv", sep=""))
      
      ## MA-plot:
      pdf(paste(WorkingDir, 
                "/Plots/", 
                GetDate(), "_", 
                tmp.name,
                ".pdf",
                sep = ""
      ), paper = "a4r")
      par(bty="n", mfrow=c(1,1))
      smoothScatter_TopTags(TT = current.out, MAIN = tmp.name)
      dev.off()
      
      #########################################################################################################
      
      rm(current.results)
      
    }
    
    ## save environment, simple dummy file name, should alter be renamed properly to match the script.R name                                           
    if (Save.Image == TRUE){
      save.image(paste(WorkingDir, "/R/", GetDate(), "_", ImageName, ".Rdata" ,sep=""))
      
    }
    
}
