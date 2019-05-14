## Template for tximport from salmon:

library(tximport)

## tx2gene list:
TX2Gene <- fread("/Volumes/Rumpelkammer/Genomes/mm10/Gencode_M20/gencode.vM20.Tx2Gene.txt", header = T, data.table = F)

## List files:
tmp.files <- list.files(path = list.dirs(path = "~/IMTB/Fischer2019/RNA-seq/salmons", recursive=FALSE, full.names = T), 
                        pattern = "quant.sf", 
                        full.names = T)

## Reformat names:
names(tmp.files) <- sapply(strsplit(tmp.files, split = "salmons|_salmon"), function(x) x[2])

## import tx:
txi <- tximport(files = tmp.files, type = "salmon", txIn = T, txOut = F, tx2gene = TX2GENE)
