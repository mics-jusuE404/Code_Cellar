#!/bin/bash

## Call peaks with Genrich both on replicate groups indicated by BASENAME_rep*_dedup.bam

########################################################################################################################

## mm10 consensus blacklist for ATAC-seq:
Blacklist="/scratch/tmp/a_toen03/Genomes/mm10/mm10_consensusBL.bed"

## Sort BAM by name:
if [[ ! -e ./GenrichDir ]]; then
  mkdir GenrichDir; fi

########################################################################################################################

## Genrich for groups:
function GENRICHGROUP {
  
  BASENAME=$1
  MODE=$3
  
  ## Only if replicates are present:
  if [[ $(ls ${BASENAME}*dedup.bam | wc -l) < 2 ]]; then
    echo '[WARNING]: Only one sample found for' $BASENAME '-- skipping group-level peak calling'
    exit 0
    fi
  
  FILES=$(ls ${BASENAME}*dedup.bam | xargs | awk '{gsub(" ", ",");print}')
  
  if [[ ${MODE} == "PE" ]]; then
    ls ${BASENAME}*dedup.bam | parallel "samtools sort -n -@ 4 -m 1G -o ./GenrichDir/{} {}"
    cd GenrichDir
    Genrich -E $2 -t $FILES -j -l 200 -q 0.01 -o ${BASENAME}_peaks.narrowPeak
    fi
    
  if [[ ${MODE} == "SE" ]]; then  
    Genrich -S -E $2 -t $FILES -w 100 -j -l 200 -q 0.01 -o ./GenrichDir/${BASENAME}_peaks.narrowPeak
    fi
    
}; export -f GENRICHGROUP

########################################################################################################################

## Genrich for individual samples:
function GENRICHSINGLE {
  
  BASENAME=$1
  MODE=$3
  
  if [[ ${MODE} == "PE" ]]; then
    ls ${BASENAME}*dedup.bam | parallel "samtools sort -n -@ 3 -m 1G -o ./GenrichDir/{} {}"
    cd GenrichDir
    Genrich -E $2 -t ${BASENAME}_dedup.bam -j -l 200 -q 0.05 -o ${BASENAME}_peaks.narrowPeak
    fi
    
  if [[ ${MODE} == "SE" ]]; then  
    Genrich -S -E $2 -t ${BASENAME}_dedup.bam -w 100 -j -l 200 -q 0.05 -o ./GenrichDir/${BASENAME}_peaks.narrowPeak
    fi
  
}; export -f GENRICHSINGLE

########################################################################################################################

## Example for single-sample in PE mode:
ls *_dedup.bam \
  | awk -F "_dedup.bam" '{print $1 | "sort -u"}' \
  | parallel "GENRICHSINGLE {} $Blacklist PE 2> {}_genrich.log"

