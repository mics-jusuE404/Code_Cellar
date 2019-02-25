## Load Salmon quantifications by tximport:

setwd("/scratch/tmp/a_toen03/RNA-seq/")

################################################################################################################################################

require(tximport)
require(data.table)

################################################################################################################################################

## Gencode v28 Tx2Gene
TX2Gene <- fread("/scratch/tmp/a_toen03/Genomes/hg38/Gencode_v28/gencode.v28.TX2Gene.txt", header = T, data.table = F)

## Table that indicates which SRRs belong to one patient
Patient   <- fread("path/to/dir/table1.tsv", data.table = F, header = T)
Cellline  <- fread("path/to/dir/table2.tsv", data.table = F, header = F)


################################################################################################################################################

## Tximportfiles:
TxImport_AT <- function(QuantDir){
  
  ## Get the paths for all the quantification files and then rename them to only get the SRR-ID,
  ## assuming they are all in a folder (...)/fastq and have name SRR(...)_salmonK${IDXLENGTH}:
  tmp.files <- list.files(path = list.dirs(QuantDir, recursive=FALSE), pattern = "quant.sf", full.names = T)
  names(tmp.files) <- sapply(strsplit(tmp.files, split = "/fastq/|_salmon"), function(x) x[2])
  names(tmp.files) <- gsub("/", "", names(tmp.files))
  return(tmp.files)
}

## Patients:
Salmons_Patients <- TxImport_AT(QuantDir = "./Patients/fastq/")
Salmons_Lines    <- TxImport_AT(QuantDir = "./Cell_Lines/fastq/")


################################################################################################################################################

## Import everthing:
Txi_all <- tximport(files = unlist(sapply(ls(pattern = "Salmons_"), get)),
                    type = "salmon",
                    txIn = T, txOut = F, 
                    tx2gene = TX2Gene)
