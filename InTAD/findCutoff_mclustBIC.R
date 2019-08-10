## Check for a data-driven expression cutoff in RNA-seq data.
## Code adapted from filterGeneExpr.R from InTAD:

findCutoff_mclustBIC <- function(count.matrix, do.log2 = TRUE, log2.prior = 1, plotExprDistr = FALSE){
  
  require(mclust)
  require(InTAD)
  
  ## function from InTAD::filterGeneExpr.R
  get.enr.bg.normfit <- function(x) {
    
    result <- c(NA, NA)
    
    # NOTE: no idea how not to depend on mclust
    # requireNamespace() did not work properly
    
    if (requireNamespace("mclust")) {
      names(result) <- c("mean","sd")
      
      # fit two normal distributions to data
      # model <- Mclust(na.omit(x), G=2, modelNames="V")
      model <- mclust::Mclust(na.omit(x), G=2)
      
      # bg corresponds to first (lower mean) normal
      result[1] <- model$parameters$mean[1]
      result[2] <- sqrt(model$parameters$variance$sigmasq[1])
      
      # store additional information in attributes
      attr(result, "proportions") <- model$parameters$pro
      attr(result, "enr.mean") <- as.numeric(model$parameters$mean[2])
      attr(result, "enr.sd") <- sqrt(model$parameters$variance$sigmasq[2])
      
      # return results
    } else {
      message("Mclust library not detected, skipping cut value adjustment...")
    }
    result
  }
  
  if (do.log2) count.matrix <- log2(count.matrix + log2.prior)
  
  vals <- as.vector(count.matrix)
  pcut <- 0.01
  
  d1 <- density(vals,na.rm=TRUE)
  bg1 <- get.enr.bg.normfit(vals)
  
  if (sum(is.na(bg1)) == 0) {
    ratiocut <- qnorm(pcut,  mean=bg1[1], sd=bg1[2], lower.tail=FALSE)
    
    if (plotExprDistr) {
      bgProp1 = attr(bg1, "proportions")[1]
      bgProp2 = attr(bg1, "proportions")[2]
      plot(d1, col='gray', lwd=3, main="Normal Mixture Model",
           ylim=c(0,2.5),xlab="Expression", ylab="Density")
      lines(x=d1$x,
            y=dnorm(d1$x, mean=bg1[1], sd=bg1[2])*bgProp1,
            lwd=3, col="blue")
      lines(x=d1$x,
            y=dnorm(d1$x, mean=attr(bg1, "enr.mean"),
                    sd=attr(bg1, "enr.sd")) * bgProp2,
            lwd=3, lty=2, col="green")
      abline(v=ratiocut, lwd=3, col="red")
      legend(x="topright", bty="n", lty=c(1,1,2,1,NA),
             lwd=c(3,3,3,3,NA),
             col=c("gray","blue","green","red","black"),
             legend=c(sprintf("alldata (n=%d)", d1$n),
                      sprintf("background (%.1f%%)", 100*bgProp1),
                      sprintf("foreground (%.1f%%)", 100*bgProp2),
                      sprintf("P = %.3g (n_sig=%d)",pcut,
                              length(which(vals>=ratiocut)))))
    }
    
    cutVal <- ratiocut
    
  }
  return(cutVal)
}
