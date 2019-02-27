## Load Salmon quantifications by tximport:

setwd("/scratch/tmp/a_toen03/RNA-seq/")

################################################################################################################################################

require(tximport)
require(data.table)

################################################################################################################################################

## Gencode v28 Tx2Gene
TX2Gene <- fread("/scratch/tmp/a_toen03/Genomes/hg38/Gencode_v28/gencode.v28.TX2Gene.txt", header = T, data.table = F)

## Assuming the *_salmon folders are in a folder termed <Salmon>:
TxImport_Custom <- function(QuantDir){
  
  ## List folders and get names:
  message("Checking for quantification files")
  tmp.files <- list.files(path = list.dirs(path = QuantDir, recursive=FALSE), pattern = "quant.sf", full.names = T)
  names(tmp.files) <- sapply(strsplit(tmp.files, split = "/|_salmon"), function(x) x[3])
  
  ## Do tximport:
  if ( exists("TX2Gene") == FALSE) stop("TX2Gene does not exist")
  message("Importing quantifications")
  return(
    tximport(files = tmp.files,
             type = "salmon",
             txIn = T, txOut = F, 
             tx2gene = TX2Gene)
  )
  
}

################################################################################################################################################

ImmGenRNA_tximport <- TxImport_Custom(QuantDir = "RNA-seq/Salmons")
