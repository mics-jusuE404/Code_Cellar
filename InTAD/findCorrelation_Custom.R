## This is the InTAD::findCorrelations.R script, modified to enable chunk-wise processing
## of signal/gene pairs to reduce memory consumption.

## This function is unchanged but renamed it to avoid conflicts
findGeneCorrelation_custom <- function(x, signalVals, countVals, corMethod) {
    genes <- x$geneid
    pId <- as.character(unique(x$peakid))
    
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
        corDat[[i]] <- data.frame(gene=genes[i], #eudis= euDis,
                                corr=cors$estimate,
                                pvalue=cors$p.value,
                                stringsAsFactors =FALSE)
    }
    ##corDF <- do.call(rbind, lapply(corDat, function(x) x))
    corDF <- data.table::rbindlist(corDat)
    return(corDF)
}

## Custom function to perform correlation analysis splitted over (chunk.sizes) elements per go
## to reduce memory consumption
findCorrelation_custom <- function(object,               ## InTADSig object
                                   chunk.size = 5000,    ## number of elements per chunk /default to be used on the 1.5TB node!)
                                   method = "pearson",   ## correlation method
                                   Studyname = NULL,     ## studyname for output to disk ./InTAD_raw/InTAD_raw_Studyname_permutcycle_chunk(n).tsv
                                   current.cycle = NULL, ## from the lapply function that makes the permutations the current cycle \
                                   total.cycle = NULL){  ## only needed to print correct file name

    if (is.null(Studyname)) Studyname <- "blankStudy"
    
    if (!is(object, "InTADSig"))
        stop("Object must be an InTADSig!")

    if (length(object@signalConnections)  == 0)
        stop("No signals and genes are combined in TAD!
            Use combineInTAD() function.")
    
    ## Couple of filtering steps unchanged from original:
    sigCons <- object@signalConnections
    txs <- rowRanges(object@sigMAE[["exprs"]])
    allIdGnX <- sigCons[!sapply(sigCons, is.null)]
    allIdGnY <- allIdGnX[sapply(allIdGnX, function(x) dim(x)[1]) > 0]
    assay.signal <- assay(object@sigMAE[["signals"]])
    assay.exprs  <- assay(object@sigMAE[["exprs"]])
    
    ## This function is adapted from original script and makes a handy output from allRes
    ## but uses rbindlist instead of do.call(rbind, lists.of.dataframes) to be a bit faster
    AfterCor <- function(allRes=allRes){
      
      allRes <- lapply(allRes, function(x) na.omit(x))
      
      allX <- list()
      for(i in seq_len(length(allRes))) {
        gnname <- values(txs)$gene_name[match(allRes[[i]]$gene,
                                              values(txs)$gene_id)]
        rnk <- allRes[[i]]$gene[order(allRes[[i]]$corr, decreasing=TRUE)]
        corRank <- match(allRes[[i]]$gene, rnk)
        
        peakName <- names(allRes)[i]
        tad <- unique(allIdGnX[[peakName]]$tad)
        
        allX[[i]] <- data.frame(peakid=rep(peakName, length(gnname)),
                                tad=rep(tad,length(gnname)),
                                gene=allRes[[i]]$gene, name=gnname,
                                cor=allRes[[i]]$corr,
                                pvalue=allRes[[i]]$pvalue,
                                corRank=corRank,
                                stringsAsFactors =FALSE)
      }
      return(rbindlist(allX))
    }
    
    ## launch the main function:
    total.chunks <- ceiling(length(allIdGnY) / chunk.size)
    trash.can <- lapply(X = seq(1, total.chunks), function(CHUNK){
             
     if (CHUNK == 1){
       tmp.lower <<- 1
       tmp.upper <<- chunk.size
     }
             
     message(paste("[Info]: Processing chunk",CHUNK,"of",total.chunks,"for cycle",current.cycle,"of",total.cycle))
      
     if (tmp.upper > length(allIdGnY)) tmp.upper <<- length(allIdGnY)
     allIdGnY.use <- allIdGnY[names(allIdGnY)[c(seq(tmp.lower, tmp.upper))]]
     
     allRes <- parallel::mclapply(allIdGnY.use,
                                  findGeneCorrelation_custom,
                                  assay.signal,
                                  assay.exprs,
                                  method,
                                  mc.cores=detectCores()/2)
              
      ## Initialize next round:
      tmp.lower <<- tmp.upper + 1
      tmp.upper <<- tmp.lower + chunk.size - 1
             
      ## if last chunk set tmp.upper differently
      if (CHUNK == total.chunks) tmp.upper <- allIdGnY
             
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
    
      } ## end lapply-loop
    )
    message("[Info]: Done processing all chunks")
}
