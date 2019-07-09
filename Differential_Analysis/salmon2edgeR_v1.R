## Wrapper for tximport followed by edgeR (only quality control and CPM calculation, no DEG).
## See README.md for details.

Salmon2edgeR <- function(SalmonDir,           ## the path to the folder with the samlom quantifications called <salmons>
                         GlobalName,          ## a prefix to store all variables in
                         Tx2Gene = "",        ## the tab-delim list with transcript2gene conversions
                         FilterSmallRNAs = T, ## remove from tximport all smallRNAs as these are typically not well-captured
                         FilterFile = "",     ## if "" use the internal default or specify a new file with gene names to be filtered
                         plotMAall = F        ## T/F/none to plot all possible combination, averaged per group or nothing
                         
){
  GetDate <- function(){ format(Sys.Date(), "%Y%m%d") }
  ####################################################################################################################################################
  ## Check if all packages are installed:
  tmp.check <- sapply(c("tximport", "data.table", "edgeR", "csaw"), function(x) x%in%installed.packages())
  if ( (FALSE %in% unlist(tmp.check)) == TRUE ) {
    tmp.which <- which(tmp.check == FALSE)
    stop(paste("These packages are not installed [ ", paste(names(tmp.which), collapse = " | "), " ]", sep=""))
  }
  
  ####################################################################################################################################################
  
  library(tximport)
  library(data.table)
  library(edgeR)
  library(csaw)
  
  ####################################################################################################################################################
  
  ## Use the default tx2gene list if not specified explicitely:
  if (Tx2Gene == "") {
    TX2Gene <- fread("/Volumes/Rumpelkammer/Genomes/mm10/Gencode_M20/gencode.vM20.Tx2Gene.txt", header = T, data.table = F)
  }
  
  ## Use the default filtering file list if not specified explicitely:
  if (FilterFile == "") {
    FilterFile <- fread("/Volumes/Rumpelkammer/Genomes/mm10/Gencode_M20/Filtered_Files/gene_name_smallRNA_TEC.txt", header = F, data.table = F)
  }
  
  ## List quantification files:
  tmp.files <- list.files(path = list.dirs(path = SalmonDir, recursive=FALSE, full.names = T), 
                          pattern = "quant.sf", 
                          full.names = T)
  
  ## Reformat names to only leave file name without any path:
  names(tmp.files) <- gsub("\\//|\\/", "",sapply(strsplit(tmp.files, split = "salmons|_salmon"), function(x) x[2]))
  
  ## import tx:
  txi <- tximport(files = tmp.files, type = "salmon", txIn = T, txOut = F, tx2gene = TX2Gene)
  
  ####################################################################################################################################################
  
  ## Filter smallRNAs  
  if (FilterSmallRNAs == T){
    message("Removing smallRNAs is set to TRUE")
    ## keep non-smallRNAs:
    keep <- !(rownames(txi$counts) %in% FilterFile[,1])
    tmp.names <- names(txi)
    
    txi<-lapply(txi, function(x) {
      if(class(x) == "character") return(x)
      return(x[keep,])
    })
    
  } else {
    message("Removing smallRNAs is set to FALSE")
  }
  
  ####################################################################################################################################################
  
  ## load into edgeR pretty much following the tximport manual and https://support.bioconductor.org/p/121087/ to get ength-offset corrected CPMs:
  message("Calculating offsets from tximport")
  cts <- txi$counts
  normMat <- txi$length
  normMat <- normMat/exp(rowMeans(log(normMat)))
  o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
  y <- DGEList(cts)
  ## y would now be ready for filterByExpr() followed by dispersion estimation
  y <- scaleOffset(y, t(t(log(normMat)) + o)) 

  ## export CPMs:
  message("Saving offset-corrected CPMs")
  se <- SummarizedExperiment(assays = y$counts)
  names(assays(se))[1] <- "counts"
  se$totals <- y$samples$lib.size
  assay(se, "offset") <- y$offset
  se.cpm <- calculateCPM(se, use.norm.factors = F, use.offsets = T, log = F)
  
  ## Save everything as variable:
  assign(x = paste(GlobalName, ".DGElist", sep=""), value = y, envir = .GlobalEnv)
  ## and to disk:
  if (!dir.exists("./Lists")) dir.create("./Lists")
  write.table(x = data.frame(GeneName=rownames(se.cpm), se.cpm), quote = F, row.names = F, col.names = T, sep="\t", 
              file = paste("./Lists/", paste(GetDate(), GlobalName, "CPM.tsv", sep="_"), sep = ""))
  
  ####################################################################################################################################################
  
  PLOTDIR <- "./Plots"
  if (!dir.exists("./Plots")) dir.create("./Plots")
  
  ## MDS plot:
  mds <- plotMDS(x = calculateCPM(se, use.norm.factors = F, use.offsets = T, log = T), bty="n", cex = 0.8, plot = F)
  pdf(file = paste(PLOTDIR, "/", GetDate(), "_", GlobalName,"_MDSplot.pdf", sep="") , onefile=T, paper='A4r') 
  plot(x = mds$x, y = mds$y, bty = "n", main = paste("MDS plot of ", GlobalName, sep=""), 
       xlim = c(floor(min(mds$x)), ceiling(max(mds$x))), ylim = c(floor(min(mds$y)), ceiling(max(mds$y))),
       pch = 20, cex=0.8)
  text(mds$x, mds$y, labels=names(mds$x), cex= 0.7, pos=3)
  dev.off()
  assign( paste(GlobalName, "_MDS", sep=""), mds, envir = .GlobalEnv)
  
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
  
  if (plotMAall != "none"){
    
    if (plotMAall == F){
      message("Printing MA plots averaged over replicates")
      
      if ("FALSE" %in% grepl("rep", colnames(se.cpm))) message("Warning: Not all sample names end with _rep* -- skipping MAplot part!")
      
      ## If naming convention *_rep* is ok, proceed:
      if( length( grep("FALSE", grepl("rep", colnames(se.cpm))) ) == 0){
        
        ## Average all replicates per condition:
        names.unique <- unique(sapply(strsplit(colnames(se.cpm), split="_rep"), function(x) x[1]))
        
        ## Make all pairwise comparisons:
        tmp.combn <- combn(x = names.unique, m = 2)
        
        pdf(file = paste(PLOTDIR, "/", GetDate(), "_", GlobalName,"_MAplots.pdf", sep="") , onefile=T, paper='A4') 
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
    
    if (plotMAall == T){
      
      message("Printing MA plots for all sample combinations")
      
      tmp.combn <- combn(x = colnames(se.cpm), m = 2)
      
      pdf(file = paste(PLOTDIR, "/", GetDate(), "_", GlobalName,"_MAplots.pdf", sep="") , onefile=T, paper='A4') 
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
  
  ####################################################################################################################################################
  
}

## Example when leaving everything at default:
## Salmon2edgeR(SalmonDir = "~/path/to/salmons/", plotMAall = F, GlobalName = "Testname")  
