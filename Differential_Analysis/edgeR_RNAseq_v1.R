## Template for differential analysis of RNA-seq data with edgeR assuming that salmon2edgeR.R was run before.

run_edgeR <- function(DGElist,        ## the DGElist produced by salmon2edgeR_v*.R
                      GlobalName,     ## suffix for variables
                      Coldata,        ## Coldata specifying the levels for Design
                      Design="",      ## Design parameter
                      Contrasts="",   ## Contrasts to test
                      FilterByExpr=T, ## use FilterByExpression function
                      GlmTreat.Fc=""  ## minimum fold change to test against
  ){ 
  
  GetDate <- function(){ gsub("^20", "", format(Sys.Date(), "%Y%m%d")) }
  
  ########################################################################################################################
  
  packageS <- c("edgeR")
  if (length(grep("FALSE", (packageS %in% rownames(installed.packages())))) > 0){
    stop("Package(s): ", packageS[which( packageS %in% rownames(installed.packages()) == "FALSE")], " are not installed!")
  }
  
  ########################################################################################################################
  
  library(edgeR)

  ########################################################################################################################
  
  ## define samples via coldata:
  y <- DGElist
  y$group <- Coldata
  
  ## apply recommended filter:
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
    
  ## Test all specified contrasts:
  message("Testing all Contrasts (total of ", dim(Contrasts)[2],")")
    
  for (i in seq(1,dim(Contrasts)[2])){
      
    ## current contrast:
    if (GlmTreat.Fc == ""){
      message("Null hypothesis is FC=0")
      current.results <- glmQLFTest(fit, contrast = Contrasts[,i])
    }
      
    if (GlmTreat.Fc != "" && is.numeric(GlmTreat.Fc)){
      message("Null hypothesis is FC = ",GlmTreat.Fc, " for ", attr(Contrasts, "dimnames")$Contrasts[i])
      current.results <- glmTreat(glmfit = fit, contrast = Contrasts[,i], lfc = log2(GlmTreat.Fc))
    }
      
    ## Save the FDR-adjusted TT:
    current.out <- topTags(current.results, n=Inf, adjust.method="BH", sort.by="none")
    current.out <- current.out$table
    assign( paste(GlobalName, "_topTags_", gsub("-", "_", attr(Contrasts, "dimnames")$Contrasts[i]), sep=""),
            current.out, envir = .GlobalEnv)
    rm(current.results)
  }
  
  assign(paste(GlobalName, "_final.DGElist", sep=""), y, envir = .GlobalEnv)
  assign(paste(GlobalName, ".Fit", sep=""), fit, envir = .GlobalEnv)
  
}

########################################################################################################################
########################################################################################################################

## Example:
## coldata   <- data.frame(SampleName = (...),
##                         Condition  = (...))

## design    <- model.matrix(~ 0 + Condition, data=coldata)

## contrasts <- makeContrasts(C1=(conditionB-conditionA, levels = design)

## run_edgeR(DGElist = foo.DGElist, GlobalName = "foo", 
##           Coldata = coldata, 
##           Design = design, 
##           Contrasts = contrasts)

