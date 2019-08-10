## This is the InTAD::findCorrelations.R script, modified to enable chunk-wise processing
## of signal/gene pairs to reduce memory consumption.

## This function is unchanged but renamed it to avoid conflicts
findGeneCorrelation_custom <- function(x, signalVals, countVals, corMethod) {
    genes <- x$geneid
    pId <- as.character(unique(x$peakid))
    #print(pId)

    corDat <- list()
    for (i in seq_len(length(genes))) {
        if (! (as.character(genes[i]) %in% row.names(countVals))) {
            #print(paste("NOT FOUND!", genes[i]))
            next
        }
        # log2 is taken into account in the generation of object
        countSel <- as.numeric( countVals[as.character(genes[i]), ])
        #print(as.character(genes[i]) )
        if ( sum(countSel) == 0) {
            next
        }
        sigSel = as.numeric(signalVals[pId,])
        cors <- suppressWarnings(cor.test( sigSel, countSel, method=corMethod))
        euDis = as.numeric(dist(rbind ( sigSel, countSel)))
        corDat[[i]] <- data.frame(gene=genes[i], eudis= euDis,
                                corr=cors$estimate,
                                pvalue=cors$p.value,
                                stringsAsFactors =FALSE)
    }
    corDF <- do.call(rbind, lapply(corDat, function(x) x))
    return(corDF)
}

## Custom function to perform correlation analysis splitted over chunk.sizes elements per go
## to reduce memory consumption
findCorrelation_custom <- function(object,               ## InTADSig object
                                   chunk.size = 1000,    ## number of elements per chunk   
                                   method = "pearson",   ## correlation method
                                   Studyname,            ## studyname for output to disk ./InTAD_raw/InTAD_raw_Studyname_permutcycle_chunk(n).tsv
                                   current.cycle = NULL, ## from the lapply function that makes the permutations the current cycle
                                                         ## only needed to print correct file name
                                   total.cycle = NULL,
                                   Cores = 16
                                   )
  {

    if (!is(object, "InTADSig"))
        stop("Object must be an InTADSig!")

    if (length(object@signalConnections)  == 0)
        stop("No signals and genes are combined in TAD!
            Use combineInTAD() function.")
    
    
    
    sigCons <- object@signalConnections
    txs <- rowRanges(object@sigMAE[["exprs"]])
    #this null filtering was already performed in previous step
    allIdGnX <- sigCons[!sapply(sigCons, is.null)]
    allIdGnY <- allIdGnX[sapply(allIdGnX, function(x) dim(x)[1]) > 0]
    ## </turn off (ATpoint)>
    
    ## Modification (ATpoint), always do parallel processing and chunk-wise:
    
    assay.signal <- assay(object@sigMAE[["signals"]])
    assay.exprs  <- assay(object@sigMAE[["exprs"]])
    
    ## This function does the postprocessing of the raw correlation output and is basically
    ## same as in original script but removed plotting and qvalue parts.
    ## q-value calculation will be done later manually
    AfterCor <- function(allRes=allRes){
      
      allRes <- lapply(allRes, function(x) na.omit(x))
      
      allX <- list()
      for(i in seq_len(length(allRes))) {
        gnname <- values(txs)$gene_name[match(allRes[[i]]$gene,
                                              values(txs)$gene_id)]
        rnk <- allRes[[i]]$gene[order(allRes[[i]]$corr, decreasing=TRUE)]
        corRank <- match(allRes[[i]]$gene, rnk)
        #cluster <- cltab$cluster[match(names(allRes)[i], cltab$id)]
        peakName <- names(allRes)[i]
        tad <- unique(allIdGnX[[peakName]]$tad)
        #tad <- inTADpairs$tad[match(names(allRes)[i], inTADpairs$peakid)]
        allX[[i]] <- data.frame(peakid=rep(peakName, length(gnname)),
                                tad=rep(tad,length(gnname)),
                                gene=allRes[[i]]$gene, name=gnname,
                                cor=allRes[[i]]$corr,
                                pvalue=allRes[[i]]$pvalue,
                                eucDist = allRes[[i]]$eudis,
                                corRank=corRank,
                                stringsAsFactors =FALSE)
      }
      
      return( do.call(rbind, lapply(allX, function(x) x)) )
    }
    
    ## Initialize by deviding into how many 1000er-chunks the process is split,
    ## if fewer than 1000 present set to 1:
    tmp.lower <- 1
    tmp.upper <- chunk.size
    total.chunks <- ceiling(length(allIdGnY) / tmp.upper)
    trash.can <- lapply(X = seq(1, total.chunks),
           FUN = function(CHUNK){
             
             message(
               paste("[Info]: Processing chunk", 
                     CHUNK, "of", total.chunks, 
                     "for cycle", current.cycle, "of", total.cycle)
               )
             
             if (tmp.upper > length(allIdGnY)) tmp.upper <- length(allIdGnY)
        
             allIdGnY.use <- allIdGnY[names(allIdGnY)[c(seq(tmp.lower, tmp.upper))]]
      
             allRes <- parallel::mclapply(allIdGnY.use,
                                          findGeneCorrelation_custom,
                                 assay.signal,
                                 assay.exprs,
                                 method,
                                 mc.cores=Cores)
              
             ## Initialize next round:
             tmp.lower <- tmp.upper + 1
             tmp.upper <- tmp.upper*2
      
             ## post processing
             tmp.result <- AfterCor(allRes = allRes)
             
             if (CHUNK == 1) true.false <- TRUE
             if (CHUNK > 1)  true.false <- FALSE
             
             if (is.null(current.cycle)) current.cycle <- 1
             
             ## write results to avoid keeping it in memory
             if (!dir.exists("./InTAD_raw")) dir.create("./InTAD_raw")
             write.table(x = tmp.result, 
                         append = !true.false, 
                         file = paste0("./InTAD_raw/", "InTAD_raw_", Studyname, "_permut", current.cycle, "_chunk_", CHUNK, ".tsv"),
                         quote = FALSE,
                         sep = "\t",
                         row.names = FALSE,
                         col.names = true.false)
             rm(tmp.result, allRes)
    
      } ## end for-loop
    )
    message("[Info]: Done processing all chunks")
}
