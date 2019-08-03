#!/bin/bash

## Script produces RLE-normalized bigwig files for each sample in the folder with
## RLE calculated over the a sample group where ${Basename}_foo_rep*_dedup.bam 
## everything in front of _rep is considered one group.
## Assumes that ./MacsDir contains ${Basename}_foo-stringent.bed containing reference peaks.
## Will merge the stringent peaks per group, make a count matrix, run DESeq2 for SF estimation,
## feed the SF into awk to correct the raw mosdepth bedGraphs with it.
## Will then average all RLE bigwigs between replicates.

## First draft, works ok but lacks all kinds of control steps to stop if tools or files are misisng.

Rscript="$HOME/anaconda3_things/anaconda3/envs/R_env/bin/Rscript"
Wiggle="$HOME/anaconda3_things/anaconda3/envs/wiggle/bin/wiggletools"

####################################################################################################################################

cat <<EOF > sizeFactors.R
## Script reads a count matrix from stdin (assuming featureCounts format from SAF),
## and spills out size factors from DESeq2 with $1 = name and $2 = sf to stdout:

packageS <- c("DESeq2", "data.table")
if (length(grep("FALSE", (packageS %in% rownames(installed.packages())))) > 0){
  stop("Package(s): ", packageS[which( packageS %in% rownames(installed.packages()) == "FALSE")], " are not installed!")
}

suppressMessages(library(DESeq2))
suppressMessages(library(data.table))

## read data:
counts <- fread('cat /dev/stdin', skip = 1, header = T, data.table = F)

if (ncol(counts) == 7) {
  write.table(data.frame(c("ONESAMPLE")), sep="\t", col.names = F, row.names = T,
              quote = F, file = "./sizeFactors.txt")
  stop("Only one sample present, so skipping calculation!")              
}

## get rowMeans
counts <- as.matrix(counts[,7:ncol(counts)])

## size factors:
SizeFactors <- estimateSizeFactorsForMatrix(counts = counts)

names(SizeFactors) <- gsub("_dedup.bam", "", names(SizeFactors))

## write:
write.table(data.frame(SizeFactors), sep="\t", col.names = F, row.names = T,
            quote = F, file = stdout())
EOF

####################################################################################################################################

function SF {
  
  Rscript=$2
  
  ## the basename for the overall group, e.g. PU1 for PU1_<condition>_<rep>_dedup.bam
  Basename=$1
    
  ## consensus peaks:
  ls ./MacsDir/${Basename}_*stringent.bed \
    | xargs cat | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    | tee ./MacsDir/${Basename}_merged_stringentPeaks.bed \
    | awk 'OFS="\t" {print $1"_"$2+1"_"$3, $1, $2+1, $3, "+"}' > ./MacsDir/${Basename}_merged_stringentPeaks.saf
  
  ## count matrix:  
  featureCounts \
    -a ./MacsDir/${Basename}_merged_stringentPeaks.saf \
    -F SAF \
    -T 8 \
    -o ./MacsDir/${Basename}_merged_stringentPeaks.countmatrix \
    ${Basename}*dedup.bam
    
  ## size factors:
  cat ./MacsDir/${Basename}_merged_stringentPeaks.countmatrix | "${Rscript}" sizeFactors.R
  
}; export -f SF

ls *_dedup.bam | awk -F "_rep" '{print $1 | "sort -u"}' | parallel -j 6 "SF {} ${Rscript} >> sizeFactors.txt"

####################################################################################################################################

## Simple part: Mosdepth for everything:
ls *_dedup.bam | awk -F "_dedup.bam" '{print $1}' | parallel -j 15 "mosdepth -t 4 {} {}_dedup.bam"

####################################################################################################################################

## Take the raw mosdepth output and make normalized bigwig based on the size factors:
function BigwigSF {

  Basename=$1
    
  ## bamcoverage multiplies with the size factor but deseq2 divides so invert the deseq factor:
  if [[ $(cat sizeFactors.txt | wc -l | xargs) > 1 ]]; then
    FACTOR=$(grep ${Basename} sizeFactors.txt | cut -f2)
    fi
    
  if [[ $(cat sizeFactors.txt | wc -l | xargs) == 1 ]]; then
    echo '[Error]: Only one single sample in sizeFactors.txt'
    exit 1
    fi
  
  ## Bigwig with 80bp fragments to smooth the signal:
  bgzip -c -d -@ 2 ${Basename}.per-base.bed.gz \
  | mawk -v SF=${FACTOR} 'OFS="\t" {print $1, $2, $3, $4/SF}' \
  | bg2bw -i /dev/stdin -c tmp_chromSizes.txt -o ${Basename}_rle.bigwig
  rm \
    ${Basename}.per-base* \
    ${Basename}.mosdepth*
  
}; export -f BigwigSF

ls *dedup.bam | awk -F "_dedup.bam" '{print $1}' | parallel -j 32 "BigwigSF {}"

####################################################################################################################################

## Now use wiggletools
function BigWigMean {
  
  ## basename here split at <_rep>:
  Basename=$1
  Wiggle=$2
  
  Files=$(ls ${Basename}*rep*rle.bigwig)
  
  ${Wiggle} mean ${Files} \
  | wigToBigWig stdin tmp_chromSizes.txt ${Basename}_mean_rle.bigwig
 
}; export -f BigWigMean

ls *_rep*rle.bigwig | awk -F "_rep" '{print $1}' | parallel -j 4 "BigWigMean {} ${Wiggle}"
