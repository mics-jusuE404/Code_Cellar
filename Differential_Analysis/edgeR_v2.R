## Template for differential analysis of RNA-seq data with edgeR starting from tximport files:

run_edgeR <- function(TXI,         ## tximport output
                      COLDATA,     ## coldata specifying the levels for design
                      DESIGN,      ## design parameter
                      CONTRASTS,   ## contrasts to test
                      NAME,        ## name assigned to this analysis
                      plotMAall=F, ## T/F for plotting MAs for all combinations or only aggregated per group
                      WORKINGDIR){ ## wdir for the plots to be saved
  
  ########################################################################################################################
  packageS <- c("tximport", "csaw", "statmod", "edgeR")
  if (length(grep("FALSE", (packageS %in% rownames(installed.packages())))) > 0){
    stop("Package(s): ", packageS[which( packageS %in% rownames(installed.packages()) == "FALSE")], " are not installed!")
  }
  
  ########################################################################################################################
  library(tximport)
  library(edgeR)
  library(csaw)
  options(scipen=999)
  
  message("")
  
  ########################################################################################################################
  if (dir.exists(WORKINGDIR)){
    message("Enter working directory ", WORKINGDIR)
    setwd(WORKINGDIR)
  }
  
  if (!dir.exists(WORKINGDIR)){
    message("Creating working directory ", WORKINGDIR)
    dir.create(WORKINGDIR, showWarnings = T)
    setwd(WORKINGDIR)
  }
  
  ########################################################################################################################
  ## Code to make use of the tximport information, 100% copied from tximport vignette:
  message("Calculating offsets from tximport")
  cts <- TXI$counts
  normMat <- TXI$length
  normMat <- normMat/exp(rowMeans(log(normMat)))
  
  o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
  y <- DGEList(cts)
  
  y <- scaleOffset(y, t(t(log(normMat)) + o))
  
  assign( paste(NAME, "_unfilteredDGElist", sep=""), y, envir = .GlobalEnv)
  
  # filtering
  keep <- filterByExpr(y)
  y <- y[keep, ]
  
  ## Custom part: Export normalized counts following https://support.bioconductor.org/p/121087/,
  ## making ose of the offsets from tximport to correct for transcript length and depth/composition.
  ## Thanks to Aaron Lun and Mike Love for the help:
  message("Saving offset-corrected CPMs")
  se <- SummarizedExperiment(assays = y$counts)
  names(assays(se))[1] <- "counts"
  se$totals <- y$samples$lib.size
  assay(se, "offset") <- y$offset
  se.cpm <- calculateCPM(se, use.norm.factors = F, use.offsets = T, log = F)
  
  ## Export:
  assign( paste(NAME, "_CPM", sep=""), se.cpm, envir = .GlobalEnv)
  
  ## Draw MDS:
  mds <- plotMDS(x = calculateCPM(se, use.norm.factors = F, use.offsets = T, log = T), bty="n", cex = 0.8, plot = F)
  
  pdf(file = paste(NAME,"_MDSplot.pdf", sep="") , onefile=T, paper='A4r') 
  plot(x = mds$x, y = mds$y, bty = "n", main = paste("MDS plot of ", NAME, sep=""), 
       xlim = c(floor(min(mds$x)), ceiling(max(mds$x))), ylim = c(floor(min(mds$y)), ceiling(max(mds$y))),
       pch = 20, cex=0.8)
  text(mds$x, mds$y, labels=names(mds$x), cex= 0.7, pos=3)
  dev.off()
  assign( paste(NAME, "_MDS", sep=""), mds, envir = .GlobalEnv)
  
  ########################################################################################################################
  ## Define samples via COLDATA:
  y$group <- COLDATA
  
  ## dispersion estimate and model fit:
  message("Gene-wise dispersion estimates")
  y <- estimateDisp(y, design = DESIGN)
  
  message("GLM fitting")
  fit <- glmQLFit(y, design = DESIGN)
  
  ## Test all specified contrasts:
  message("Testing all contrasts")
  
  for (i in 1:ncol(CONTRASTS)){
    
    current.FTest <- glmQLFTest(fit, contrast = CONTRASTS[,i])
    assign( paste(NAME, "_glmQLFTest_", attr(contr, "dimnames")$Contrasts[i], sep=""), current.FTest, envir = .GlobalEnv)
    
    current.TT    <- topTags(current.FTest, p.value = 1, n = Inf, sort.by = "none")
    assign( paste(NAME, "_topTags_", attr(contr, "dimnames")$Contrasts[i], sep=""), current.TT, envir = .GlobalEnv)
  }
  
  ########################################################################################################################
  ## Produce MA plots, using the average per replicate group:
  plotMA_custom <- function(COUNTS, MAIN = ""){
    
    R=COUNTS[,1]
    G=COUNTS[,2]  
    
    ## get mean of counts and log2FC
    M <- log2(R/G)
    A <- 0.5*log2(R*G)
    
    ## Decide y-axis limits based on rounded quantiles
    YLIM <- c(floor(quantile(M, 0.0001, na.rm = T)), ceiling(quantile(M, 0.9999, na.rm = T)))
    
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
  
  if (plotMAall == "N"){
    if (plotMAall == F){
      message("Printing MA plots averaged over replicates")
      
      if ("FALSE" %in% grepl("rep", colnames(se.cpm))) message("Warning: Not all sample names end with _rep* -- skipping MAplot part!")
      
      ## If naming convention *_rep* is ok, proceed:
      if( length( grep("FALSE", grepl("rep", colnames(se.cpm))) ) == 0){
        
        ## Average all replicates per condition:
        names.unique <- unique(sapply(strsplit(colnames(se.cpm), split="_rep"), function(x) x[1]))
        
        ## Make all pairwise comparisons:
        tmp.combn <- combn(x = names.unique, m = 2)
        
        pdf(file = paste(NAME,"_MAplots.pdf", sep="") , onefile=T, paper='A4') 
        par(mfrow=c(1,1), bty="n")
        for (i in 1:ncol(tmp.combn)){
          
          one <- se.cpm[,grep(tmp.combn[1,i], colnames(se.cpm))]
          two <- se.cpm[,grep(tmp.combn[2,i], colnames(se.cpm))]
          
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
      
      tmp.combn <- combn(x = colnames(se.cpm), m = 2)
      
      pdf(file = paste(NAME,"_MAplots.pdf", sep="") , onefile=T, paper='A4') 
      par(mfrow=c(1,1), bty="n")
      for (i in 1:ncol(tmp.combn)){
        
        one <- se.cpm[,grep(tmp.combn[1,i], colnames(se.cpm))]
        two <- se.cpm[,grep(tmp.combn[2,i], colnames(se.cpm))]
        
        ## MA:
        tmp.df <- data.frame(one, two)
        colnames(tmp.df) <- tmp.combn[,i]
        
        
        plotMA_custom(COUNTS = tmp.df, MAIN = paste("MA-plot: ", paste(tmp.combn[,i], collapse=" vs "), sep=""))
        
      }; suppressMessages(dev.off())
    }
  }
}


########################################################################################################################
########################################################################################################################

## Example:
#coldata <- data.frame(NAME = colnames(txi$counts),
#                      GENOTYPE = c(rep("KO", 5), rep("WT", 6)),
#                      STIMULUS = c(rep("Ca", 2), rep("U", 3), rep("Ca", 3), rep("U", 3)),
#                      FACTORIZED = c(rep("A", 2), rep("B", 3), rep("C", 3), rep("D", 3)))

#design <- model.matrix(~ 0 + FACTORIZED, data=coldata)

#contr <- makeContrasts( KOuKOCa = FACTORIZEDB-FACTORIZEDA,
#                       WUCaWTu = FACTORIZEDD-FACTORIZEDC ,
#                        levels = design)

#run_edgeR(TXI = txi, COLDATA = coldata, DESIGN = design, CONTRASTS = contr, NAME = "test", plotMAall = "N", WORKINGDIR = "~/MUH")

