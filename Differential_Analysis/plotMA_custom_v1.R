## Simple wrapper script to make MA-plots from counts in a matrix or dataframe.
## Usage example at the bottom of the script.
## Automatically chooses y-axis limits based on rounded quantiles unless user overwrites this with YLIM option.
## Can use smoothScatter or heatScatter for plotting. The latter requires the LSD package.
## smoothScatter is faster and produces smaller (in terms of disk space) images,
## heatscatter is somewhat more fancy but produces quite large images that do not really fit into Powerpoint.

plotMA_custom <- function(COUNTS,           ## 2-column data.frame or matrix with the counts to plot
                          LOGGED = F,       ## F=counts are not logged, T=counts are already log2
                          MAIN = "",        ## plot title
                          MODE = "smooth",  ## smooth for smoothScatter plot or heatscatter for LSD::heatscatter
                          YLIM="",          ## user-defined ylims, if empty will be chosen based on quantiles
                          REPLACE.ZERO = 1){ ## uses the specified number (default 1) to replace zeros to avoid log2(0) 
  
par(bty="n")
## Check if LSD is installed:
if (MODE == "heatscatter" && !("LSD" %in% rownames(installed.packages()))){
  message("LSD package is not installed, using smoothScatter instead")
  MODE <- "smooth"
}

## F if provided counts are on the regular (=non-log) scale
if (LOGGED == F){
  R=COUNTS[,1]
  G=COUNTS[,2]  
  
  R[which(R == 0)] <- REPLACE.ZERO
  G[which(G == 0)] <- REPLACE.ZERO
}

## T if counts are already on log2 scale
if (LOGGED == T){
  R=2^COUNTS[,1]
  G=2^COUNTS[,2]  
}

## get mean of counts and log2FC
M <- log2(R/G)
A <- 0.5*log2(R*G)

## Decide y-axis limits based on rounded quantiles
if (YLIM[1] == ""){
  YLIM <- c(floor(quantile(M, 0.0001, na.rm=T)), ceiling(quantile(M, 0.9999, na.rm=T)))
}

## using smoothScatter
if (MODE == "smooth"){
  smoothScatter(A, M, main = MAIN, bty="n",
                xlab="mean of normalized counts",
                ylab="log2FC", ylim=YLIM)
  abline(h=0)  
}

## using heatscatter (more fancy but pdfs will be much larger)
if (MODE == "heatscatter"){
  library(LSD)
  heatscatter(A, M, main = MAIN, bty="n", pch=20,
              xlab="mean of normalized counts",
              ylab="log2FC", ylim=YLIM)
  abline(h=0)  
}

## plot points beyond y-axis limits as triangles to the y-axis limits
if (MODE == "smooth") COL = "darkblue"
if (MODE == "heatscatter") COL = "firebrick"
points(x = A[which(M < YLIM[1])], y = rep(YLIM[1], length(which(M < YLIM[1]))), pch=17, cex = 0.7, col=COL)
points(x = A[which(M > YLIM[2])], y = rep(YLIM[2], length(which(M > YLIM[2]))), pch=17, cex = 0.7, col=COL)

}

## Example usage:
par(mfrow=c(2,1))
plotMA_custom(COUNTS = cpm1[,c(1:2)], LOGGED = F, MAIN = "default y-axis limits based on quantiles", MODE = "smooth")
plotMA_custom(COUNTS = cpm1[,c(1:2)], LOGGED = F, MAIN = "user-defined y-axis (-1,1)", MODE = "smooth", YLIM = c(-1,1))
