## Basic volcano plot, autoadjusting x- and y-limits

Volcano_AT <- function(DFresults, FC="logFC", pval="FDR", 
                       colBase = "grey85", colDown = "darkgoldenrod", colUp = "darkblue", 
                       main = "Volcano", p.thres = 0.05, fc.thres = 0){
  
  ## x/y dimensions for plot:
  y.max <- 2* round ( round( -log10( min(DFresults[,which(colnames(DFresults)==pval)]) ) ) / 2)
  x.max <- ceiling(max(DFresults[,which(colnames(DFresults)==FC)]))
  x.min <- floor(min(DFresults[,which(colnames(DFresults)==FC)]))
  
  ## add singificant up:
  which.up <- which( DFresults[,which(colnames(DFresults)==pval)] < p.thres &
                       DFresults[,which(colnames(DFresults)==FC)] > fc.thres )
  which.dw <- which( DFresults[,which(colnames(DFresults)==pval)] < p.thres &
                       DFresults[,which(colnames(DFresults)==FC)] < fc.thres )
  
  ## Main plot in grey:
  plot(x = DFresults[,which(colnames(DFresults)==FC)],
       y = -log10( DFresults[,which(colnames(DFresults)==pval)] ),
       pch = 20, col = colBase,
       xlim=c(x.min, x.max),
       ylim=c(0, y.max),
       main = paste(main, " -- [", "Down: ", length(which.dw), "] -- ", "[Up: ", length(which.up), "]", sep=""),
       cex = 0.6,
       xlab="log2FC", ylab=paste("-log10(",pval,")", sep=""), bty="n")
  
  points(y = -log10( DFresults[which.up, which(colnames(DFresults)==pval)]),
         x = DFresults[which.up, which(colnames(DFresults)==FC)],
         col=colUp, pch=20, cex=0.6)
  
  points(y = -log10( DFresults[which.dw, which(colnames(DFresults)==pval)]),
         x = DFresults[which.dw, which(colnames(DFresults)==FC)],
         col=colDown, pch=20, cex=0.6)
}
