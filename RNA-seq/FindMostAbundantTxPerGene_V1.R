## Naive approach to find most abundant transcript per gene for one or more samples/replicates,
## based on TPMs from salmon transcript quantifications. 
## Returns a data.frame with the Tx name for each replicate group or sample.

FindMostAbundantTxPerGene <- function(SalmonDir,             ## directory where salmon output folders are, expected to follow *_rep*_salmon
                                      AggregateReps = TRUE,  ## get most abundant tx per gene based on replicate average TPM
                                      Tx2Gene = NULL){       ## if T use default mm10 gencode v20
  
  #########################################################################################################################################
  
  library(data.table)
  
  #########################################################################################################################################
  
  ## Default mm10 Gencode v20, matching ensembl Tx and gene IDs:
  if (is.null(Tx2Gene)) {
    ensmus_tx2gene <- fread("/Volumes/Rumpelkammer/Genomes/mm10/Gencode_M20/ENSMUST_2_ENSMUSG.txt", 
                            data.table = F, 
                            header = F)
  }
  
  #########################################################################################################################################
  
  ## Scan salmon dirs for replicates:
  if (AggregateReps == TRUE){
    
    tmp.names <- list.files(SalmonDir, pattern = "_salmon$")  
    unique.basenames <- as.character(unique(sapply(strsplit(tmp.names, split="_rep"), function(x)x[1])))
    
    #######################################################################################################################################
    
    tmp.mostAbundant <- sapply(unique.basenames, function(x){
      
      ## list all files that have the same basename:
      tmp.reps <- grep(paste0("^", x, "_rep"), list.files(SalmonDir, pattern = "_salmon$"), value = TRUE)
      
      ## get the tx-gene-tpm table for each replicate of a group:
      tmp.byReplicates <- 
      lapply(tmp.reps, function(m){
        tmp.read  <- fread(paste0(SalmonDir, "/", m, "/quant.sf"), data.table = F, header = T, select = c(1,4))
        tmp.merge <- merge(ensmus_tx2gene, tmp.read, by.x = "V1", by.y = "Name")
      })
      
      tmp.df <- do.call(cbind, tmp.byReplicates)
      tmp.df <- data.frame(tmp.df[,c(1,2, grep('TPM', colnames(tmp.df)))])
      
      if ( length(grep('TPM', colnames(tmp.df))) > 1){
        tmp.df <- data.frame(tmp.df[c(1,2)], TPM = rowMeans(tmp.df[,seq(3, ncol(tmp.df))]))  
      }
      
      ## keep most abundant (by TPM) transcript per gene:
      tmp.df <- tmp.df[order(-tmp.df$TPM) , ] #sort by id and reverse of abs(value)
      return( tmp.df[ !duplicated(tmp.df[,2]), 1] ) 
      
    })
    
    
    #######################################################################################################################################
    
  }
  
  return( data.frame( tmp.mostAbundant ))
  
}
  
a<- FindMostAbundantTxPerGene(SalmonDir = "~/IMTB/Project_Leukemia/RNAseq_July2019/salmons/")
  
