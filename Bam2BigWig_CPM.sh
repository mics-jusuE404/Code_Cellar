#!/bin/bash

######################################################################################################
##
## Take an indexed BAM file and create a CPM-normalized bigwig , using:
## - mosdepth for a compressed bedgraph (.bed.gz),
## - bgzip for decompression
## - mawk for custom to-CPM conversion and
## - bg2bw (git: cancerit/cgpBigWig) to read the bg from stdin and output the bw
##
######################################################################################################

BAM=$1

######################################################################################################

command -v mosdepth >/dev/null 2>&1 || { echo >&2 "[ERROR]: mosdepth is not in PATH"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "[ERROR]: samtools is not in PATH"; exit 1; }
command -v mawk >/dev/null 2>&1 || { echo >&2 "[ERROR]: mawk is not in PATH"; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo >&2 "[ERROR]: bgzip is not in PATH"; exit 1; }
command -v bc >/dev/null 2>&1 || { echo >&2 "[ERROR]: bgzip is not in PATH"; exit 1; }

######################################################################################################

## Check bam file:

if [[ $# -ne 1 ]] ; then
    echo 'Usage: ./Bam2BigWig_CPM.sh in.bam'
    exit 1
fi

DO_EXIT=$(echo 'BAM file looks corrputed -- exiting' && exit 1)
samtools quickcheck $1 && echo '[INFO]: Processing' $1 || ${DO_EXIT}

## Check if indexed:
if [[ ! -f ${1}.bai ]]; then
  echo '[INFO]: BAM file is not indexed -- indexing now:'
  sambamba index -t 8 $1
  fi

######################################################################################################

## Functions:
SCALE_FACTOR=$(bc <<< "scale=8;1000000/$(samtools idxstats $1 | awk '{SUM+=$3} END {print SUM}')") 

## Depth with mosdepth:
mosdepth ${1%.bam} $1

## ChromSizes:
samtools idxstats $1 | \
  mawk 'OFS="\t" {print $1, $2 | "sort -k1,1 -k2,2n"}' | \
    grep -v '*' > ${1%.bam}_chromsizes.txt

## Normalize CPM:
bgzip -c -d ${1%.bam}.per-base.bed.gz | \
  mawk -v SF=${SCALE_FACTOR} 'OFS="\t" {print $1, $2, $3, $4*SF}' | \
  bg2bw -i /dev/stdin \
        -c ${1%.bam}_chromsizes.txt \
        -o ${1%.bam}_CPM.bigwig  
    
######################################################################################################

rm ${1%.bam}.per-base*
rm ${1%.bam}_chromsizes.txt
rm ${1%.bam}.mosdepth*
