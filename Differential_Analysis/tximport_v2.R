## Function to import quantifications from salmon stored in a folder "salmons"
## By default will remove smallRNA genes as they should not be captured well during library
## prep and therefore eventually only (probably) inflate multiple testing burden.
## If everything is at default, use tx2gene and filtering file from my local machine for gencode mouse:

Run_Tximport <- function(Filepath,            ## the path to the folder with the samlom quantifications called <salmons>
                         Tx2Gene = "",        ## the tab-delim list with transcript2gene conversions
                         FilterSmallRNAs = T, ## remove from tximport all smallRNAs as these are typically not well-captured
                         FilterFile = "")     ## if "" use the internal default or specify a new file with gene names to be filtered
                         {
  
  library(tximport)
  library(data.table)
  
  ## Use the default mouse gencode tx2gene if not specified:
  if (Tx2Gene == "") {
    TX2Gene <- fread("/Volumes/Rumpelkammer/Genomes/mm10/Gencode_M20/gencode.vM20.Tx2Gene.txt", header = T, data.table = F)
  }
  
  if (FilterFile == "") {
    FilterFile <- fread("/Volumes/Rumpelkammer/Genomes/mm10/Gencode_M20/Filtered_Files/gene_name_smallRNA_TEC.txt", header = F, data.table = F)
  }
  
  ## List files ti import:
  tmp.files <- list.files(path = list.dirs(path = Filepath, recursive=FALSE, full.names = T), 
                          pattern = "quant.sf", 
                          full.names = T)
  
  ## Reformat names to only leave file name without any path:
  names(tmp.files) <- gsub("\\//|\\/", "",sapply(strsplit(tmp.files, split = "salmons|_salmon"), function(x) x[2]))
  
  ## import tx:
  txi <- tximport(files = tmp.files, type = "salmon", txIn = T, txOut = F, tx2gene = TX2Gene)

  ## Filter smallRNAs  
  if (FilterSmallRNAs == T){
    message("Removing smallRNAs is set to TRUE")
    ## keep non-smallRNAs:
    keep <- !(rownames(txi$counts) %in% FilterFile[,1])
    tmp.names <- names(txi)
    
    txi_filtered<-lapply(txi, function(x) {
      if(class(x) == "character") return(x)
      return(x[keep,])
    })
    return(txi_filtered)
  } else {
    message("Removing smallRNAs is set to FALSE")
    return(txi)
  }
      
}

## Example when leaving everything at default:
Run_Tximport(FILEPATH = "~/path/to/salmons/")  
