## Wrapper for differential analysis with edgeR downstream of timport from salmon quantifications.
## Expected naming convention is that the quantification folder (where quant.sf is in) suffixes with _salmon
GetDate <- function(){ gsub("^20", "", format(Sys.Date(), "%Y%m%d")) }

salmon2edgeR <- function(Basename,                   ## Prefix for all elements saved to disk/envir
                         SalmonDir,                  ## the path to the folder with the salmon quantifications
                         WorkingDir = "./",          ## default working directory
                         Tx2Gene,                    ## the tab-delim list with transcript2gene conversions
                         FilterSmallRNAs = TRUE ,    ## remove from tximport all smallRNAs as these are typically not well-captured
                         SmallRNAFile,               ## the file that contains the genes to filter out
                         Save.RawDGElist = FALSE,    ## if TRUE saves the DGElist prior to all FilterByExpr etc
                         MakeMAplots = "groupwise",  ## one of c("all", "groupwise", "none") for explaroatory MA-plots
                         MakePCAplot = TRUE,         ## guess what
                         Return.PCA = FALSE,         ## whether to return the PCA data
                         log.PCA = TRUE,             ## PCA based on log2(+1) CPM
                         Return.SE = FALSE,          ## whether to return the summarized experiments with norm. counts used for PCA
                         nPCAgenes = 1000,           ## the number of genes to select for PCA
                         Return.tximport = F,        ## whether to save the tximport result to envir
                         Coldata = NULL,             ## Coldata for differential expression
                         Design = NULL,              ## design for edgeR
                         Contrasts = NULL,           ## contrasts based on makeContrasts
                         Thresh.glmTreat = NULL,     ## fold change to test against as null hypothesis     
                         FilterByExpr = TRUE,        ## whether or not to use that filter from edgeR
                         Save.ImageName = NULL       ## if not null, save image to that file in the format ./R/GetDate()_<name>.Rdata
                         
){
  
  GetDate <- function(){ gsub("^20", "", format(Sys.Date(), "%Y%m%d")) }
  
  ####################################################################################################################################################
  
  ## Check if all packages are installed:
  tmp.check <- sapply(c("tximport", "data.table", "edgeR", "csaw", "statmod", "DESeq2"), function(x) x%in%installed.packages())
  if ( (FALSE %in% unlist(tmp.check)) == TRUE ) {
    tmp.which <- which(tmp.check == FALSE)
    stop(paste("These packages are not installed [ ", paste(names(tmp.which), collapse = " | "), " ]", sep=""))
  }
  
  ####################################################################################################################################################
  
  suppressPackageStartupMessages(require(tximport))
  suppressPackageStartupMessages(require(data.table))
  suppressPackageStartupMessages(require(edgeR))
  suppressPackageStartupMessages(require(csaw))
  suppressPackageStartupMessages(require(statmod))
  suppressPackageStartupMessages(require(DESeq2))
  
  ####################################################################################################################################################
  
  smoothScatter_TopTags <- function(TT, SIG.THRESHOLD=0.05, YLIM="", MAIN=NULL) {
    
    ## ylim based on quantiles:
    if (YLIM[1] == ""){
      YLIM <- c(floor(quantile(TT$logFC, 0.0001, na.rm=T)), ceiling(quantile(TT$logFC, 0.9999, na.rm=T)))
    }
    
    par(bty="n")
    if (is.null(MAIN)) t.main <- deparse(substitute(TT))
    if (!is.null(MAIN)) t.main <- MAIN
    
    smoothScatter(x = TT$logCPM, y = TT$logFC,
                  xlab="average logCPM", 
                  ylab="log2FC", 
                  main=paste(t.main,
                             " at ", 
                             "FDR = ", 
                             SIG.THRESHOLD*100, 
                             "%", sep=""), 
                  ylim=YLIM, cex.main=0.75)
    
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
  
  
  ####################################################################################################################################################
  
  ## Create the specified working directory:
  if (!dir.exists(WorkingDir)) dir.create(WorkingDir, recursive = T)
  setwd(WorkingDir)
  
  if (!dir.exists("./Plots/")) dir.create("./Plots/")
  PLOTDIR <- "./Plots/"
  
  ####################################################################################################################################################
  
  ## Use the default tx2gene (mm10, gencode) list if not specified explicitely:
  TX2Gene <- fread(Tx2Gene, header = F, data.table = F, sep = "\t")
  
  SmallRNAFile <- fread(SmallRNAFile, header = F, data.table = F, sep = "\t")
  
  
  ## List quantification files:
  tmp.files <- list.files(path = list.dirs(path = SalmonDir, recursive=FALSE, full.names = T), 
                          pattern = "quant.sf", 
                          full.names = T)
  
  ## Check if quantification directories end with "_salmon":
  tmp.check <- sapply(strsplit(sapply(strsplit(tmp.files, split="/"), function(x)rev(x)[2]), split="_"), function(x)rev(x)[1])
  if (length(tmp.check) != sum(grepl("salmon", tmp.check))){
    warning(paste(paste0("Not all quantification directories in ", SalmonDir, " appear to end with _salmon."),
                  paste0("That might mess up correct colname formatting. Better check directory names!"), sep="\n"))
  }
  
  ## Remove the directory name and the mandatory suffix of the quantification directory:
  names(tmp.files) <- gsub("\\//|\\/", "",sapply(strsplit(tmp.files, 
                                                          split = paste(rev(strsplit(SalmonDir, split="/")[[1]])[1], "_salmon", sep="|")), 
                                                 function(x) x[2]))
  
  ####################################################################################################################################################
  
  ## Aggregate transcript level estimates to the gene level using tximport:
  txi <- tximport(files = tmp.files, type = "salmon", txIn = T, txOut = F, tx2gene = TX2Gene)
  
  ####################################################################################################################################################
  
  ## Filter smallRNAs as in standard RNA-seq they probably represent artifacts and are not reliable (size selection 200bp during lib. prep):
  if (FilterSmallRNAs == T){
    message("Removing smallRNAs is set to TRUE")
    ## keep non-smallRNAs:
    keep <- !(rownames(txi$counts) %in% SmallRNAFile[,1])
    tmp.names <- names(txi)
    
    txi<-lapply(txi, function(x) {
      if(class(x) == "character") return(x)
      return(x[keep,])
    })
    
  } else {
    message("Removing smallRNAs is set to FALSE")
  }
  
  if (Return.tximport == TRUE) assign(paste(Basename, ".tximport", sep=""), txi, envir = .GlobalEnv)
  
  ####################################################################################################################################################
  
  ## Calculate normalized counts corrected with the offsets from tximport for length, 
  ## see => https://support.bioconductor.org/p/121087/ for details.
  message("Calculating offsets from tximport")
  cts <- txi$counts
  normMat <- txi$length
  normMat <- normMat/exp(rowMeans(log(normMat)))
  o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
  y <- DGEList(cts)
  ## y would now be ready for filterByExpr() followed by dispersion estimation
  y <- scaleOffset(y, t(t(log(normMat)) + o)) 
  
  ## if T save the DGElist prior to any further filtering:
  if (Save.RawDGElist == TRUE) assign(paste0(Basename, "_raw.DGElist"), y, envir = .GlobalEnv)
  
  ## export CPMs:
  message("Saving offset-corrected CPMs")
  se <- SummarizedExperiment(assays = y$counts)
  names(assays(se))[1] <- "counts"
  se$totals <- y$samples$lib.size
  assay(se, "offset") <- y$offset
  se.cpm <- calculateCPM(se, use.norm.factors = FALSE, use.offsets = TRUE, log = FALSE)
  
  ## To envir and disk:
  assign(x = paste(Basename, "_CPM.df", sep=""), value = data.frame(se.cpm), envir = .GlobalEnv)
  
  if (!dir.exists("./Lists")) dir.create("./Lists")
  write.table(x = data.frame(Gene=rownames(se.cpm), se.cpm), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t", 
              file = paste("./Lists/", paste(GetDate(), Basename, "CPM.tsv", sep="_"), sep = ""))
  
  ####################################################################################################################################################
  
  if (!dir.exists("./Plots")) dir.create("./Plots")
  
  ## PCA plot using DESeq2 template:
  if (MakePCAplot == TRUE){
    
    message("Plotting PCA")
    
    ## New SE object for later use with DESeq2::plotPCA:
    se.deseq2 <- SummarizedExperiment(assays = se.cpm)
    se.deseq2$samples <- as.factor(colnames(assay(se.deseq2)))
    se.deseq2$groups  <- as.factor(sapply(strsplit(colnames(assay(se.deseq2)), split="_rep"), function(x)x[1]))
    
    if(log.PCA) assay(se.deseq2) <- log2(assay(se.deseq2)+1)
    
    if (Return.SE) assign(paste0(Basename, "_cpm.se"), se.deseq2, envir = .GlobalEnv)
    
    se.deseq2 <- DESeqTransform(se.deseq2)
    
    pdf(paper = "a4", file = paste0(WorkingDir, "/Plots/", GetDate(), "_", Basename, "_PCA_top", round(nPCAgenes/1000, digits=2), "k.pdf"))
    print(DESeq2::plotPCA(object = se.deseq2, intgroup="groups", nPCAgenes))
    dev.off()
    
    if (Return.PCA) {
      assign(paste0(Basename, ".PCAdata"),
             DESeq2::plotPCA(object = se.deseq2, intgroup="groups", nPCAgenes, returnData = TRUE),
             envir = .GlobalEnv)
    }
    
  }
  
  ####################################################################################################################################################
  
  ## MA plots to get an idea of the data:
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
  
  ## get normalized counts
  
  if (MakeMAplots != "none"){
    
    if (MakeMAplots == "groupwise"){
      message("Printing MA plots averaged over replicates")
      
      if ("FALSE" %in% grepl("rep", colnames(se.cpm))) message("Warning: Not all sample names end with _rep* -- skipping MAplot part!")
      
      ## If naming convention *_rep* is ok, proceed:
      if( length( grep("FALSE", grepl("rep", colnames(se.cpm))) ) == 0){
        
        ## Average all replicates per condition:
        names.unique <- unique(sapply(strsplit(colnames(se.cpm), split="_rep"), function(x) x[1]))
        
        ## Make all pairwise comparisons:
        tmp.combn <- combn(x = names.unique, m = 2)
        
        pdf(file = paste(PLOTDIR, "/", GetDate(), "_", Basename,"_MAplots.pdf", sep="") , onefile=T, paper='A4') 
        par(mfrow=c(1,1), bty="n")
        for (i in 1:ncol(tmp.combn)){
          
          one <- se.cpm[,grep(tmp.combn[1,i], colnames(se.cpm))]
          two <- se.cpm[,grep(tmp.combn[2,i], colnames(se.cpm))]
          
          if (ncol(data.frame(one)) == 1) average.one <- one
          if (ncol(data.frame(one)) > 1) average.one <- rowMeans(one)
          if (ncol(data.frame(two)) == 1) average.two <- two
          if (ncol(data.frame(two)) > 1) average.two <- rowMeans(two)
          
          ## MA:
          tmp.df <- data.frame(average.one, average.two)
          colnames(tmp.df) <- tmp.combn[,i]
          
          
          plotMA_custom(COUNTS = tmp.df, MAIN = paste("Sample: ", paste(tmp.combn[,i], collapse=" vs "), sep=""))
          
        }; dev.off()
      }
    }
    
    if (MakeMAplots == "all"){
      
      message("Printing MA plots for all sample combinations")
      
      tmp.combn <- combn(x = colnames(se.cpm), m = 2)
      
      pdf(file = paste(PLOTDIR, "/", GetDate(), "_", Basename,"_MAplots.pdf", sep="") , onefile=T, paper='A4') 
      par(mfrow=c(1,1), bty="n")
      for (i in 1:ncol(tmp.combn)){
        
        one <- se.cpm[,grep(tmp.combn[1,i], colnames(se.cpm))]
        two <- se.cpm[,grep(tmp.combn[2,i], colnames(se.cpm))]
        
        ## MA:
        tmp.df <- data.frame(one, two)
        colnames(tmp.df) <- tmp.combn[,i]
        
        
        plotMA_custom(COUNTS = tmp.df, MAIN = paste(paste(tmp.combn[,i], collapse=" over "), sep=""))
        
      }; suppressMessages(dev.off())
    }
  }
  
  ####################################################################################################################################################
  
  #### Start edgeR part, stop here if no coldata/design is given:
  
  if (is.null(Coldata) | is.null(Design)) {
    message("No Design given therefore stop before differential analysis part")
    message("[Ended]:", "\t", Basename)
    #https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
    stop_quietly <- function() {
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      stop()
    }; stop_quietly()
  }
  
  message("Differential analysis")
  
  ## add coldata to DGElist:
  y$group <- Coldata
  
  ## Apply filter (might make sense if one expects few DGE):
  if (FilterByExpr == T){
    keep <- filterByExpr(y, design = Design)
    y <- y[keep, , keep.lib.sizes=FALSE]
  }
  
  ## dispersion estimates:
  message("Gene-wise dispersion estimates")
  y <- estimateDisp(y, design = Design)
  
  ## GLM fitting:
  message("GLM fitting")
  fit <- glmQLFit(y, design = Design)
  assign(paste0(Basename, ".glmQLFit"), fit, envir = .GlobalEnv)
  
  ## Test all specified contrasts:
  message("Testing all Contrasts (total of ", dim(Contrasts)[2],")")
  
  for (i in seq(1,dim(Contrasts)[2])){
    
    ## current contrast:
    if (is.null(Thresh.glmTreat)){
      message("Null hypothesis is FC=0")
      current.results <- glmQLFTest(fit, contrast = Contrasts[,i])
    }
    
    if (!is.null(Thresh.glmTreat) && is.numeric(Thresh.glmTreat)){
      message("Null hypothesis is FC = ",Thresh.glmTreat, " for ", attr(Contrasts, "dimnames")$Contrasts[i])
      current.results <- glmTreat(glmfit = fit, contrast = Contrasts[,i], lfc = log2(Thresh.glmTreat))
    }
    
    ## Save the FDR-adjusted TT:
    current.out <- topTags(current.results, n=Inf, adjust.method="BH", sort.by="none")
    current.out <- current.out$table
    
    tmp.assignname <- paste0(Basename, "_topTags_", gsub("-", "_", attr(Contrasts, "dimnames")$Contrasts[i]))
    
    ## MA-plots:
    pdf(paste0("./Plots/", GetDate(), "_", tmp.assignname, ".pdf"), paper = "a4r")
    par(mfrow=c(1,1), bty="n")
    smoothScatter_TopTags(current.out, MAIN = tmp.assignname)
    dev.off()
    
    current.out <- data.frame(Gene = rownames(current.out), current.out)
    rownames(current.out) <- NULL
    
    write.table(sep="\t", quote = F, row.names = F, col.names = T,
                file = paste0("./Lists/", GetDate(), "_", tmp.assignname, ".tsv"), x = current.out)
    
    assign( tmp.assignname,
            current.out, envir = .GlobalEnv)
    
    rm(current.results)
  }
  
  assign(paste0(Basename, "_final.DGElist"), y, envir = .GlobalEnv)
  
  if (!is.null(Save.ImageName)){
    if (!dir.exists("./R")) dir.create("./R")
    save.image(paste0("./R/", GetDate(), "_", Save.ImageName, ".Rdata"))
  }

}

####################################################################################################################################################


