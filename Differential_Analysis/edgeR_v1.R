## Template for differential analysis of RNA-seq data starting from tximport files:

library(tximport)
library(edgeR)

run_edgeR <- function(TXI = txi, COLDATA = coldata, DESIGN= design, CONTRASTS = contrasts, NAME=""){
  
  if (NAME == "") stop("Please enter a sample name!")
  
  ## - TXI is the output of tximport()
  TXI <- txi
  ########################################################################################################################
  ## Correction for gene/isoform length and depth using TMM,
  ## 100% copied from tximport vignette:
  message("Library normalization")
  cts <- TXI$counts
  normMat <- TXI$length
  normMat <- normMat/exp(rowMeans(log(normMat)))

  o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
  y <- DGEList(cts)
  y <- scaleOffset(y, t(t(log(normMat)) + o))
  # filtering
  keep <- filterByExpr(y)
  y <- y[keep, ]
  # y is now ready for estimate dispersion functions see edgeR User's Guide
  
  ## Export:
  assign( paste(NAME, "_data_DGEList", sep=""), y, envir = .GlobalEnv)
  
  ########################################################################################################################
  ## Define samples via coldata:
  y$group <- coldata
  
  ## dispersion estimate and model fit:
  message("Gene-wise dispersion estimates and model fit")
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  
  ## test all contrasts:
  message("Testing all contrasts")
  
  for (i in 1:ncol(CONTRASTS)){
    
    current.FTest <- glmQLFTest(fit, contrast = CONTRASTS[,i])
    assign( paste(NAME, "_data_FTest_", attr(contr, "dimnames")$Contrasts[i], sep=""), current.FTest, envir = .GlobalEnv)
    
    current.TT    <- topTags(current.FTest, p.value = 1, n = Inf, sort.by = "none")
    assign( paste(NAME, "_data_topTags_", attr(contr, "dimnames")$Contrasts[i], sep=""), current.TT, envir = .GlobalEnv)
  }
}

########################################################################################################################
########################################################################################################################

## Example:
coldata <- data.frame(NAME = colnames(txi$counts),
                      GENOTYPE = c(rep("KO", 5), rep("WT", 6)),
                      STIMULUS = c(rep("Ca", 2), rep("U", 3), rep("Ca", 3), rep("U", 3)),
                      FACTORIZED = c(rep("A", 2), rep("B", 3), rep("C", 3), rep("D", 3)))

design <- model.matrix(~ 0 + FACTORIZED, data=coldata)

contr <- makeContrasts( KOuKOCa = FACTORIZEDB-FACTORIZEDA,
                        WUCaWTu = FACTORIZEDD-FACTORIZEDC ,
                        levels = design)

run_edgeR(TXI = txi, COLDATA = coldata, DESIGN = design, CONTRASTS = contr, NAME = "TEST")


  
