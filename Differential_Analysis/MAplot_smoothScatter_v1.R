## MAplots from edgeR's topTags or DESeq2's results objects coloring significant genes.

MAplot_smoothScatter <- function(Input,                      ## a topTags (edgeR) or results (DESeq2) object
                                 Preset = "DESeq2",          ## DESeq2 or edgeR
                                 Signif.Thresh = 0.05,       ## color genes below this padj/FDR in Signif.Color
                                 Signif.Color = "firebrick", ## guess what this option does
                                 Ylim = "",                  ## ylim, e.g. c(-4,4), if empty choose automatically
                                 Plot.Legend = T             ## whether to print a legend with summary stats
                                 )
  {         
  
  ###############################################################################################################
  
  ## Get fold change, average expression and padj/FDR from Input:
  if (Preset == "DESeq2"){
    logfc  <- Input$log2FoldChange
    logcpm <- log2(Input$baseMean+1) ## as DESeq2 does reports non-logged counts
    fdr    <- Input$padj
  }
  
  if (Preset == "edgeR"){
    logfc  <- Input$logFC
    logcpm <- Input$logCPM
    fdr    <- Input$FDR
  }
  
  ###############################################################################################################
  
  ## If no ylim given, choose based on quantiles to capture most data points:
  if (Ylim[1] == ""){
    Ylim <- c(floor(quantile(logfc, 0.0001, na.rm=T)), ceiling(quantile(logfc, 0.9999, na.rm=T)))
  }
  
  par(bty="n") ## because boxes around plots are ugly :)
  smoothScatter(x = logcpm, y = logfc,
                xlab="average expression", 
                ylab="log2FC", 
                main=paste(deparse(substitute(Input)),
                           " at ", 
                           "FDR = ", 
                           Signif.Thresh*100, 
                           "%", sep=""), 
                ylim=Ylim)
  
  ###############################################################################################################
  
  ## Color significant genes:
  signif.logfc     <- c(logfc[  which( logfc[which(fdr < Signif.Thresh)] > 0) ],
                        logfc[  which( logfc[which(fdr < Signif.Thresh)] < 0) ])
  
  signif.logcpm    <- c(logcpm[ which( logfc[which(fdr < Signif.Thresh)] > 0) ],
                        logcpm[ which( logfc[which(fdr < Signif.Thresh)] < 0) ])
  
  signif.fdr       <- c(fdr[    which( logfc[which(fdr < Signif.Thresh)] > 0) ],
                        fdr[    which( logfc[which(fdr < Signif.Thresh)] < 0) ])

  ## add significants regions in firebrick:
  points(x = signif.logcpm, y = signif.logfc, pch=20, col = Signif.Color, cex=0.1)
  
  ## add triangles for points beyond Ylim:
  lower.limit <- par("usr")[3]
  
  points(x = logcpm[which(logfc < lower.limit)], 
         y = rep( (lower.limit - lower.limit*0.025), length(which(logfc < lower.limit))),
         pch=17, cex = 0.7, col=Signif.ColorCOL) 
  
  upper.limit <- par("usr")[4]
  points(x = logcpm[which(logfc > upper.limit)], 
         y = rep( (upper.limit - upper.limit*0.025), length(which(logfc > upper.limit))),
         pch=17, cex = 0.7, col=Signif.ColorCOL) 
  
  ###############################################################################################################
  
  ## Plot legend with number of signif. UP/DOWN genes and the minimal observed FC among these genes:
  if (Plot.Legend == TRUE){
    
    legend("topright", paste( paste("UP=",sum(signif.logfc>0), sep=""), 
                              paste("DOWN=", sum(signif.logfc<0), sep=""), 
                              paste("minFC=", round(2^min(abs(signif.logfc)), digits=1), sep=""),
                              sep="\n"), bty="n") 
  }
  
  ###############################################################################################################

}
