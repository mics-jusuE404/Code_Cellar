#!/bin/bash

## First draft for ChIP-seq peak calling using replicates with Genrich,
## missing any argparse at this stage, later use flags for average length, fdr etc.

## Right now calls peaks on replicated samples at 1% FDR against mm10 blacklist without input

## mm10 consensus blacklist for ATAC-seq:
Blacklist="/scratch/tmp/a_toen03/Genomes/mm10/mm10_consensusBL.bed"

## Sort BAM by name:
if [[ ! -e ./GenrichDir ]]; then
  mkdir GenrichDir; fi

## Genrich for groups:
function GENRICH {
  
  Basename=$1
  Blacklist=$2
  Length=$3
  Group=$4
  
  ## Only if replicates are present:
  if [[ $(ls ${Basename}*dedup.bam | wc -l) < 2 ]]; then
    echo '[WARNING]: Only one sample found for' $Basename '-- skipping group-level peak calling'
    exit 0
    fi
  
  ## trigger grouped mode:
  if [[ ${Group} == "group" ]]; then
    FILES=$(ls ${Basename}*dedup.bam | xargs | awk '{gsub(" ", ",");print}')
    FDR=0.01
    fi
  
  ## trigger single mode:
  if [[ ${Group} == "single " ]]; then
    FILES=${Basename}_dedup.bam
    FDR=0.05
    fi
  
  Genrich -t $FILES -E ${Blacklist} -q ${FDR} -w ${Length} -S -o - \
  | tee ./GenrichDir/${Basename}_peaks.narrowPeak \
  | awk 'OFS="\t" {print $1"_"$2+1"_"$3, $1, $2+1, $3, "+"}' \
  > ./GenrichDir/${Basename}_peaks.saf
    
}; export -f GENRICH

ls *_dedup.bam | awk -F "_rep" '{print $1 | "sort -u"}' | parallel "GENRICH {} $Blacklist 150 group"
ls *_dedup.bam | awk -F "_dedup" '{print $1 | "sort -u"}' | parallel "GENRICH {} $Blacklist 150 single"

