## Wrapper to annotate genomic regions with annotatr and plotting of a nice pie :)

################################################################################################################################################

annotatR_Pie <- function(Ranges, 
                         Colors = NULL, 
                         Annotations = NULL, 
                         Return.Annot = TRUE, 
                         Plot.Pie = TRUE,
                         Main.Plot = NULL){
  
  
  ##############################################################################################################################################
  
  if (class(Ranges) != "GRanges") stop("Ranges must be GRanges")
  
  ##############################################################################################################################################
  
  ## Check if required packages are installed:
  packageS <- c("annotatr", "org.Mm.eg.db")
  
  if (length(grep("FALSE", (packageS %in% rownames(installed.packages())))) > 0){
    stop(call. = FALSE, "==> Package(s): ", 
         packageS[which( packageS %in% rownames(installed.packages()) == "FALSE")], 
         " are not installed!")
  }
  
  library(annotatr)
  library(org.Mm.eg.db)
  
  ##############################################################################################################################################
  
  
  ## Default mm10 basic annotations:
  
  if (is.null(Annotations)) {
    library(org.Mm.eg.db)
    tmp.annot <- suppressWarnings(build_annotations(genome = "mm10", annotations = c("mm10_genes_intergenic", "mm10_basicgenes")))
  }
  
  tmp.summarized <- data.frame(
    summarize_annotations(annotate_regions(regions = Ranges, 
                                           annotations = tmp.annot, 
                                           ignore.strand = TRUE, 
                                           quiet = TRUE))
  )
  
  if (isTRUE(Return.Annot)) return(tmp.summarized)
  
  ## Plot
  if (Plot.Pie == TRUE){

    message("Plotting pie")
    par(oma=c(0,0,0,0), mar=c(0,0,0,0))
    tmp.sum  <- sum(tmp.summarized$n)
    tmp.perc <- round( (tmp.summarized$n / tmp.sum) * 100, digits = 2)
  
    if(is.null(Colors)) {
      tmp.col <- c("firebrick", "skyblue4", "black", "darkorange3", "goldenrod3", "darkblue", "darkolivegreen4")
    }
    
    if (is.null(Main.Plot))  tmp.main  <- deparse(substitute(Ranges))
    if (!is.null(Main.Plot)) tmp.main  <- Main.Plot
    
    print(pie(tmp.perc, labels = "", col = tmp.col, lty = 1, lwd = 1, radius = 0.25))
    
    text(x=0, y = 0.34, labels = paste(tmp.main, paste0( "(", length(Ranges), " regions", ")" ) ), font = 2)
    
    tmp.text <- sapply(strsplit(tmp.summarized$annot.type, split="_"), function(x)rev(x)[1])
    tmp.text <- gsub("UTR", "'UTR", tmp.text)
    
    legend(x = -0.32, y = -0.25, 
           legend = tmp.text, fill = tmp.col, bty="n", horiz = F, 
           y.intersp = 1.4, x.intersp = 0.5, text.font=1, cex=1, ncol = 2, text.width=0.25)
    }
  
}

################################################################################################################################################
