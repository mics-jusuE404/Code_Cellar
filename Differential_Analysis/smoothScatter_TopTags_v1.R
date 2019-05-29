## Take a topTags() output from csaw/edgeR and a significance threshold and 
## produce an MA-plot with signif. regions colored in firebrick:

smoothScatter_TopTags <- function(TT, SIG.THRESHOLD=0.05, YLIM="") {
  
  ## ylim based on quantiles:
  if (YLIM[1] == ""){
    YLIM <- c(floor(quantile(TT$logFC, 0.0001, na.rm=T)), ceiling(quantile(TT$logFC, 0.9999, na.rm=T)))
  }
  
  par(bty="n")
  smoothScatter(x = TT$logCPM, y = TT$logFC,
                xlab="average logCPM", 
                ylab="log2FC", 
                main=paste(deparse(substitute(TT)),
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
