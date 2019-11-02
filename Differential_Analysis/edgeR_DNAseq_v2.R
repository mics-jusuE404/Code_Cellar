## Wrapper for edgeR/csaw-based differential analysis for DNA-seq such as ATAC/ChIP-seq.
## Minimum required input is a ranged count matrix (GRanges), an edgeR design and contrasts.
##
## First step is count normalization. For this TMM from edgeR is used based on the peak counts,
## unless Bins.gr is provided, so counts across large (e.g.10kb) bins across the genome.
## See csaw manual for details.
## Next will make PCA and MA-plots for quality- and normalization control / visualization.
## If Design = NULL, will perform no differential analysis.
## Else, will perform standard DE analysis based on the provided Design and Contrasts (makeContrasts output).
## DE tables (topTags) will be written to disk plus will make MA-plots with significant regions colored in red
##
## Depends on csaw, edgeR, DESeq2,statmod packages
##
###################################################################################################################

edgeR_DNAseq <- 
  function(Basename,                       ## the prefix for this analysis
           Counts.gr,                      ## GRanges object with raw counts
           Bins.gr = NULL,                 ## optional GRanges with counts of large bins for normalization
           Colnames.Suffix = NULL,         ## suffix to trim away from the colnames of the count matrices
           MakeMAplots = "groupwise",      ## c("all", "groupwise", "none") for MA-plots after normalization
           MakePCAplot = TRUE,             ## whether to make PCA plot based on DESeq2::plotPCA
           nTopPCA = 5000,                 ## number of top-variable elements to calculate PCs for PCA
           logPCA = TRUE,                  ## whether to use log2-transformed norm. counts for PCA
           returnDST = FALSE,              ## whether to return the DESeqTransform used for PCA as Basename.dst
           WriteCpmToDisk = TRUE,          ##
           WorkingDir = "./",              ## the directory to write files to disk
           Design = NULL,                  ## experimental design for diff. analysis
           FiltByExpr = FALSE,             ## whether or not to use that filter
           Contrasts = NULL,               ## contrasts to test
           Thresh.glmTreat = NULL,         ## minFC to test against
           Sig.Threshold = 0.05            ## FDR threshold to color significant regions in MA-plots
  )
  {
    
    GetDate <- function(){ gsub("^20", "", format(Sys.Date(), "%Y%m%d")) }
    
    ###############################################################################################################
    
    ## Check if required packages are installed:
    packageS <- c("csaw", "edgeR", "statmod", "DESeq2")
    
    if (length(grep("FALSE", (packageS %in% rownames(installed.packages())))) > 0){
      stop(call. = FALSE, "==> Package(s): ", 
           packageS[which( packageS %in% rownames(installed.packages()) == "FALSE")], 
           " are not installed!")
    }
    
    require(csaw)
    require(edgeR)
    require(statmod) ## don't ask me what this is but edgeR needs it for something
    require(DESeq2)  ## only for plotPCA
    
    ###############################################################################################################
    
    #### Some sanity checks:
    
    ## correct data format, contains count data and only numeric:
    
    CheckCountsGR <- function(GRs){
      
      if (class(GRs) != "GRanges") stop(call. = FALSE, "Counts.gr must be a GRanges object")
      if (ncol(elementMetadata(GRs)) == 0) {
        stop(call. = FALSE, "Counts.gr does not seem to consist any count columns")
      }
      if (!is.numeric(as.matrix(elementMetadata(GRs)))) {
        stop(call. = FALSE, "make sure Counts.gr countains numeric values")
      }
      
      ## Check if all samples contain the _rep information
      if(FALSE %in% unlist(lapply(colnames(elementMetadata(GRs)), function(x)grepl("_rep",x)))){
        stop(call. = FALSE, "Make sure all sample names contain _rep information, e.g. Sample1_wildtype_rep1.bam")
      }
      
    }
    
    ## Sanity check for peak counts and bibs:
    CheckCountsGR(Counts.gr)
    if (!is.null(Bins.gr)) CheckCountsGR(Counts.gr)
    
    ###############################################################################################################
    
    ## Check if MakeMAplots is one of the three supported options:
    if(!MakeMAplots %in% c("all", "groupwise", "none")){
      stop(call. = FALSE, "MakeMAplots must be one of (all, groupwise, none)")  
    }
    
    ## if MakePCAplot is T/F
    if (!is.logical(MakePCAplot) || !is.logical(logPCA)){
      warning("==> MakePCAplot and logPCA must be boolean. Settings MakePCAplot to FALSE")
      MakePCAplot <- FALSE
    }
    
    ## If not existing, create directories for plots and lists in specified working directory:
    if (!dir.exists(WorkingDir)){
      message("Creating WorkingDir: -- ", normalizePath(WorkingDir))
    }
    for (i in c("Lists", "Plots")){
      
      if (!dir.exists(paste(WorkingDir, i,sep=""))){
        dir.create(paste(WorkingDir, i,sep=""), showWarnings = TRUE)
      }
    }; rm(i)
    
    if(!is.null(Design)){
      if (!is.null(Thresh.glmTreat) & !is.numeric(Thresh.glmTreat)) stop(call. = FALSE,
                                                                         "Thresh.glmTreat is not numeric")
      if(!is.logical(FiltByExpr)) stop(call. = FALSE, "FiltByExpr must be TRUE/FALSE")
      if(is.null(Contrasts)) stop(call. = FALSE, "Please give Contrasts")
      
    }
    
    if(!is.integer(nTopPCA)) stop(call. = FALSE, "nTopPCA must be an integer")
    
    ###############################################################################################################
    
    #### Start the actual work:
    
    message("[Started]:","\t", Basename)
    
    ## Optionally remove suffix from colnames of count matrix:
    if(!is.null(Colnames.Suffix)){
      colnames(elementMetadata(Counts.gr)) <- gsub(Colnames.Suffix,"",colnames(elementMetadata(Counts.gr)))
    }
    
    ## Bring data in a format that csaw/edgeR can handle it (SummarizedExperiments):
    data.se <- SummarizedExperiment(
      list(counts=as.matrix(elementMetadata(Counts.gr))))
    
    ## add total depth (=colSums) and ranges:
    data.se$totals     <- colSums(as.matrix(elementMetadata(Counts.gr)))
    rowRanges(data.se) <- Counts.gr[,0]
    
    ###############################################################################################################
    
    #### Normalize data, either based on the peak counts or based on the large bins:
    
    ## Peaks (if no bins.gr are provided):
    if (is.null(Bins.gr)){
      message("Normalizing data based on peak counts")
      data.se <- normFactors(data.se, se.out=data.se)
    }
    
    ## Else, use bins:
    if (!is.null(Bins.gr)){
      
      message("Normalizing data based on large bins")
      
      colnames(elementMetadata(Bins.gr)) <- gsub(Colnames.Suffix, 
                                                 "", 
                                                 colnames(elementMetadata(Bins.gr)))
      
      bins.se <- SummarizedExperiment(
        list(counts=as.matrix(elementMetadata(Bins.gr))))
      
      bins.se$totals <- colSums(as.matrix(elementMetadata(Bins.gr)))
      
      nF <- normFactors(bins.se, se.out=FALSE)
      
      ## Write the norm. factors into the original data.se (= counts in peaks) object
      ## so that the bin-derived scaling factors will be used to correct the peak-based per-million counts
      data.se$norm.factors <- nF
      
    }
    
    ###############################################################################################################
    
    #### Calculate CPMs and write to disk:
    
    message("Calculating CPMs")
    CPM.Counts.gr <- Counts.gr[,0]
    
    elementMetadata(CPM.Counts.gr) <- as.matrix(
      calculateCPM(object = data.se, use.norm.factors = TRUE, log = FALSE))
    
    ## Save cpm to disk:
    if (WriteCpmToDisk == TRUE){
      message("Writing CPMs to disk")
      write.table(x = data.frame(CPM.Counts.gr)[,c(1,2,3)], 
                  file = paste(WorkingDir, "Lists/", 
                               GetDate(),"_", Basename, "_normCounts.tsv",
                               sep=""),
                  quote = FALSE, 
                  col.names = TRUE, 
                  row.names = FALSE,
                  sep="\t")
    }
    
    ###############################################################################################################
    
    #### Function to produce MA-plots with smoothScatter:
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
    
    ###############################################################################################################
    
    #### Wrapper for above function:
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
            
            plotMA_custom(COUNTS = tmp.df, MAIN = paste("Sample: ", 
                                                        paste(tmp.combn[,i], collapse=" over "), 
                                                        sep=""))
            
          }
        }
      }
      
      if (MakeMAplots == "all"){
        
        message("Printing MA-plots for all sample combinations")
        
        ## use combn to make all pairwise combinations of the samples:
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
    
    ## Call MA-plot function and append a suffix to the final pdf that indicates which norm. method was used:
    if (!is.null(Bins.gr)) SUFFIX <- "largebins"  
    if (is.null(Bins.gr))  SUFFIX <- "peakbased"
    
    if (MakeMAplots != "none"){
      
      par(mfrow=c(1,1), bty="n")
      pdf(paper = "a4", file = paste(WorkingDir,"./Plots/", GetDate(), "_MAplots_", SUFFIX, ".pdf", sep=""))
      
      Do_MAplot(CPMs = data.frame(elementMetadata(CPM.Counts.gr)))
      
      dev.off()
      
    }
    
    #########################################################################################################
    #### Plot a PCA (for this transform cpm to rangedsummarizedexperiment and add coldata):
    
    if (MakePCAplot == TRUE){
      
      message("Plotting PCA")
      
      tmp.counts <- as.matrix(elementMetadata(CPM.Counts.gr))
      
      ## by default use log2 normalized counts for PCA:
      if(logPCA) tmp.counts <- log2(tmp.counts+1)
      
      ## make RangedSummarizedExperiments object to be compatible with DESeq2::plotPCA
      rse <- SummarizedExperiment(
        list(counts=tmp.counts))
      
      rowRanges(rse) <- CPM.Counts.gr[,0]
      rse <- DESeqTransform(rse)
      
      rownames(colData(rse)) <- colnames(rse)
      colData(rse)$group <- as.factor(sapply(strsplit(colnames(rse), split="_rep"), function(x)x[1]))
      
      pdf(paper = "a4", file = paste(WorkingDir, "/Plots/", GetDate(), "_PCA_top",nTopPCA,".pdf", sep=""))
      print(plotPCA(object = rse, intgroup="group", nTopPCA))
      dev.off()
      
      if(returnDST) assign(paste0(Basename,".dst"), rse, envir = .GlobalEnv)
    }
    
    #########################################################################################################
    #### If Design is NULL stop here without differential analysis
    
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
    #### Start differential analysis:
    
    ## Bring to DGEList format:
    data.dge <- asDGEList(data.se)
    
    ## optionally filter for regions with sufficient counts for DE analysis
    if (FiltByExpr) {
      keep <- filterByExpr(data.dge, design = Design)
      #data.se <- data.se[keep,]
      data.dge <- data.dge[keep, , keep.lib.sizes=FALSE] #looks odd but is copied from csaw
    }
    
    ## dispersion estimates:
    message("Estimating dispersion")
    data.dge <- estimateDisp(data.dge, Design)
    
    ## glm data.fit:
    message("Fitting model")
    data.fit <- glmQLFit(data.dge, design = Design, robust = TRUE) #robust=T is suggested for DNA-seq in csaw
    
    ## save in environment
    assign(paste(Basename, ".DGElist", sep=""),  data.dge, envir = .GlobalEnv)
    assign(paste(Basename, ".glmQLFit", sep=""), data.fit, envir = .GlobalEnv)
    
    #########################################################################################################
    
    ## before testing of contrasts starts define a function to plot MA-plots with signif regions colored:
    smoothScatter_TopTags <- function(TT, SIG.THRESHOLD=Sig.Threshold, MAIN) {
      
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
    
    ## Test all specified contrasts, optionally against a certain FC with glmTreat,
    ## and correct for multiple testing using default BH:
    
    for (i in seq(1,dim(Contrasts)[2])){
      
      ## If Thresh.glmTreat is not specified test against null of 0:
      if (is.null(Thresh.glmTreat)){
        message("Null hypothesis is FC = 0")
        current.results <- glmQLFTest(data.fit, contrast = Contrasts[,i])
      }
      
      ## or against null of Thresh.glmTreat:
      if (!is.null(Thresh.glmTreat)){
        
        message("Null hypothesis is FC = ",
                Thresh.glmTreat, 
                " for ", 
                attr(Contrasts, "dimnames")$Contrasts[i])
        
        current.results <- glmTreat(glmfit = data.fit, 
                                    contrast = Contrasts[,i], 
                                    lfc = log2(Thresh.glmTreat))
        
      }
      
      ## Save the FDR-adjusted TopTags without sorting so input = output order:
      current.out <- topTags(current.results, n=Inf, adjust.method="BH", sort.by="none")
      current.out <- data.frame( data.frame(rowRanges(data.se))[,1:3], current.out$table)
      
      if (is.null(Thresh.glmTreat))  thr <- 0
      if (!is.null(Thresh.glmTreat)) thr <- Thresh.glmTreat
      
      ## save as variable to envir indicating the null hypothesis in the name:
      tmp.name <- paste(Basename, 
                        "_topTags_", 
                        gsub("-", "_", attr(Contrasts, "dimnames")$Contrasts[i]),
                        "_null_",
                        thr,
                        sep="")
      
      ## as variable to envir:
      assign( tmp.name,
              current.out,
              envir = .GlobalEnv)
      
      ## write results to disk
      write.table(current.out, sep="\t", col.names = T, row.names = F, quote = F,
                  file=paste(WorkingDir, "/Lists/", GetDate(), "_", tmp.name, "_NullHypo_",thr,".tsv", sep=""))
      
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
    
}
