## Up-to-date corde to read alevin output into R from the first author
## https://github.com/mikelove/tximport/issues/24

ReadAlevin <- function( base.path = NULL ){
  if (! dir.exists(base.path )){
    stop("Directory provided does not exist")
  }
  
  barcode.loc <- paste0( base.path, "alevin/quants_mat_rows.txt" )
  gene.loc <- paste0( base.path, "alevin/quants_mat_cols.txt" )
  matrix.loc <- paste0( base.path, "alevin/quants_mat.gz" )
  if (!file.exists( barcode.loc )){
    stop("Barcode file missing")
  }
  if (! file.exists(gene.loc) ){
    stop("Gene name file missing")
  }
  if (! file.exists(matrix.loc )){
    stop("Expression matrix file missing")
  }
  matrix <- as.matrix(read.csv( matrix.loc, header=FALSE))
  matrix <- t(matrix[,1:ncol(matrix)-1])
  
  cell.names <- readLines( barcode.loc )
  gene.names <- readLines( gene.loc )
  
  colnames(matrix) <- cell.names
  rownames(matrix) <- gene.names
  matrix[is.na(matrix)] <- 0
  return(matrix)
}
