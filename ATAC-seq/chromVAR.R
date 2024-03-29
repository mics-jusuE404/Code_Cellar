#!/home/a/a_toen03/anaconda3/bin/Rscript

## Basic wrapper script to run chromVAR on two-condition bulk ATAC-seq data, e.g. wildtype vs. knockout condition.
## In line 58 one specifies the peak files and on line 61 one specified the BAM files to be used for count matrix creation.
## 
##############################################################################################################################

## Packages:
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(chromVARmotifs)
library(BSgenome.Mmusculus.UCSC.mm10)
library(data.table)
library(JASPAR2018)

set.seed(2018)

## Register cores for parallel processing:
register(MulticoreParam(16))

##############################################################################################################################

## Function to load BED files as GRanges(), e.g. Bed2GR(BED = "path/to/foo.bed")
Bed2GR <- function(BED){
  tmp.bed <- fread(BED, header = F, data.table = F)
  tmp_granges <- GRanges(seqnames=as.character(tmp.bed[,1]),
                         ranges=IRanges(
                           start = tmp.bed[,2]+1, end = tmp.bed[,3])
  ); return(tmp_granges)
}

##############################################################################################################################

## Function (borrowed from chromVAR/TFBStools) to access JASPAR motifs (here the CORE vertebrate collection):
jaspar <- function (collection = "CORE", ...) 
{
  opts <- list()
  opts["tax_group"] <- "vertebrates"
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
    names(out) <- paste(names(out), TFBSTools::name(out), 
                        sep = "_")
  return(out)
}
motifs_JASPAR2018 <- jaspar()

##############################################################################################################################

## Step 1 -- Read in a list of reference ATAC-seq regions.
## The function will take narrowPeak (e.g. from macs2) and extend the summit by 200bp.
## We use 200bp to avoid excessive motif matches by chance.
peaks_200bp <- readNarrowpeaks(filename = "peaks.narrowPeak", 
                               width = 200, 
                               non_overlapping = TRUE)

## With this peak set, make a count matrix, here indicate the BAM files:
bam.list <- c("file_condition1_rep1.bam", 
              "file_condition1_rep2.bam",
              "file_condition2_rep1.bam",
              "file_condition2_rep1.bam")

## count reads per peak. colData indicates the conditions that the BAM files correspond to:
peaks_counts          <- getCounts(alignment_files = c(bam.list), 
                                   peaks           = peaks_untreated_200bp, 
                                   paired          = TRUE, 
                                   by_rg           = FALSE, 
                                   format          = "bam", 
                                   colData         = DataFrame(condition = c("condition1", 
                                                                             "condition1",
                                                                             "condition2",
                                                                             "condition2")))

##############################################################################################################################

## Step 2 -- Function that runs the standard workflow from chromVAR, which is:
## -preprocessing of the count matrix (adding GC bias, filtering for depth, all default options),
## -matching motifs to the regions
## -computeDeviations and computeVariability (output is Outname_deviation/variability)
## -the direction of the variability output is determined by checking if the deviation score 
##  is smaller or larger than zero upon the second condition.
##  If smaller, than the variability score gets a *(-1)

Bam2Deviation <- function(FragmentCounts, Peaks, Outname, Genome, motifs = motifs_JASPAR2018){
  
  ## this is pretty much copied from the chromVAR manual page:
  fragment_counts <- addGCBias(FragmentCounts, genome = Genome)
  filtered_counts <- filterSamples(fragment_counts, shiny = F)
  filtered_counts <- filterPeaks(filtered_counts)
  
  ## match motifs to the input regions.
  ## As this method does (AFAIK) not correct for multi. testing, it might be clever to use FIMO
  ## externally and feed the motif positions in manually.
  motif_ix <- matchMotifs(motifs, filtered_counts,
                          genome = Genome)
  
  ## get deviation score:
  tmp.dev <- computeDeviations(object = filtered_counts, annotations = motif_ix)
  
  assign(paste(Outname, "_deviation", sep=""),
         tmp.dev,
         envir = .GlobalEnv
  )
  
  tmp.vari <- computeVariability(object = tmp.dev)
  
  ## If the deviation from the mean is negative, so less accessable in condition1, put a minus in front of the variablility score:
  tmp.z <- assays(tmp.dev)$z
  tmp.means <- sapply(unique(tmp.dev$celltype), function(x) rowMeans(tmp.z[,grep(x, tmp.dev$celltype)]))
  tmp.direction <- as.numeric(which( tmp.means[,1] > tmp.means[,2] ))
  
  tmp.vari$variability[tmp.direction] <- tmp.vari$variability[tmp.direction] * (-1)
  
  assign(paste(Outname, "_variability", sep=""),
         tmp.vari,
         envir = .GlobalEnv
  )
}

## Run it:
Bam2Deviation(FragmentCounts = peaks_counts, Peaks = peaks_200bp, Genome = BSgenome.Mmusculus.UCSC.mm10, Outname = "FOO")                    
