## Take salmon quantifications for specified samples, average TPMs and then select the most abudant Tx per gene:

FindMostAbundantTxPerGene <- function(SalmonDir,             ## directory where salmon output folders are, expected to follow *_rep*_salmon
                                      SampleName,            ## the SampleName of the samples to be considered
                                      Tx2Gene){              ## file with $1 the ensembl tx name and $2 the ensembl gene ID without header
  
  #########################################################################################################################################
  
  library(data.table)
  
  colnames(Tx2Gene) <- c("V1", "V2")
  Tx2Gene$V1 <- sapply(strsplit(as.character(Tx2Gene$V1), split="\\."), function(x)x[1])
  Tx2Gene$V2 <- sapply(strsplit(as.character(Tx2Gene$V2), split="\\."), function(x)x[1])
  
  #########################################################################################################################################
  
  ## Scan salmon dirs for the intended data:
  tmp.names <- grep(paste0("^", SampleName), list.files(SalmonDir, pattern = "_salmon$"), value = T)  
  
  if (length(tmp.names) > 1){
    tmp.df<-do.call(cbind,
            lapply(tmp.names, function(m){
      return(fread(paste0(SalmonDir, "/", m, "/quant.sf"), data.table = F, header = T, select = c(1,4)))
    }))
    tmp.df <- tmp.df[,c(1,grep("TPM", colnames(tmp.df)))]
    tmp.df <- data.frame( ENSMUST = sapply(strsplit(as.character(tmp.df$Name), split="\\."), function(x)x[1]),
                          AverageTPM = rowMeans(tmp.df[,grep("TPM", colnames(tmp.df))]))
    
  }
  
  ## add genes:
  tmp.df <- merge(x = tmp.df, y = Tx2Gene, by.x="ENSMUST", by.y = "V1", sort = FALSE)
  
  ## keep only most abundant transcript per gene
  tmp.df <- tmp.df[order(-tmp.df$AverageTPM) , ]
  
  tmp.df <- tmp.df[ !duplicated(tmp.df$V2), ]
  
  tmp.df <- tmp.df[,c(1,3,2)]
  
  colnames(tmp.df) <- c("ENSMUST", "ENSMUSG", "TPM")
  return(tmp.df) 
  }
