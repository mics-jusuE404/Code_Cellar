#!/bin/bash

## Wrapper script to perform differential footprinting analysis on ATAC-seq data with HINT-ATAC
## (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1642-2)
## When installing from scrath, also install ghostscript and texlive core via conda.

############################################################################################################

## Paths to tools if not in $PATH
## (not in $PATH here because tool is python2.7 and I do not want any trouble with 
##  python3 and python2 things in PATH)

conda activate RGT

RGTMAIN="/home/a/a_toen03/.local/bin"
ORGANISM="mm10"
BAM1=
BAM2=

############################################################################################################

## Step-1: Do footprinting on the two groups
function Hint_Footprinting {

  BAM=$1
  PEAKS=$2
  BASENAME=$3
  ORGANISM=$4
  
  (>&2 paste -d " " <(echo '[INFO]' 'Hint-ATAC Footprinting for' "${3}" 'started on') <(date))
  
  "${RGTMAIN}"/rgt-hint footprinting \
    --atac-seq \
    --paired-end \
    --organism="${ORGANISM}" \
    --output-location=./ \
    --output-prefix="${BASENAME}" \
    "${BAM}" \
    "${PEAKS}"
    
  (>&2 paste -d " " <(echo '[INFO]' 'Hint-ATAC Footprinting for' "${3}" 'ended on') <(date))    
    
  }; export -f Hint_Footprinting

############################################################################################################

## Step-2: Motif matching

function Hint_MatchMotifs {

  BED1=$1
  BED2=$2
  BASENAME=$3
  
  (>&2 paste -d " " <(echo '[INFO]' 'Hint-ATAC MotifMatching for' "${3}" 'started on') <(date))
  
  "${RGTMAIN}"/rgt-motifanalysis matching --organism= --input-files ./match/(...) ./match/(...)
  
  (>&2 paste -d " " <(echo '[INFO]' 'Hint-ATAC MotifMatching for' "${3}" 'ended on') <(date))
  
}; export -f Hint_MatchMotifs

############################################################################################################

## Step-3: Differential Footprinting Analysis
## Don't use more than 16 cores or will complain about problems with forking whatever...
function Hint_Differential {
  
  BASENAME=$3
  
  (>&2 paste -d " " <(echo '[INFO]' 'Hint-ATAC Differential for' "${3}" 'started on') <(date))
  
  "${RGTMAIN}"/rgt-hint differential \
    --organism=mm9 \
    --bc --nc 16 \
    --mpbs-file1=./match/cDC1_mpbs.bed \
    --mpbs-file2=./match/pDC_mpbs.bed \
    --reads-file1=cDC1.bam \
    --reads-file2=pDC.bam \
    --condition1=cDC1 \
    --condition2=pDC \
    --output-location=cDC1_pDC
    
    (>&2 paste -d " " <(echo '[INFO]' 'Hint-ATAC Differential for' "${3}" 'started on') <(date))
    
}; export -f Hint_Differential 
  
############################################################################################################  
  
