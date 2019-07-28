#!/bin/bash

## Wrapper for peak calling on ChIP-seq data using PePr on replicate groups and single files:

if [[ ! -e ./PepperDir ]]; then mkdir PepperDir; fi

function PEPPER {
  
  BASENAME_IP=$1
  BASENAME_IGG=$2
  BAMTYPE=$3
  THRESH=$4
  PEAKTYPE=$5
  
  ######################################################################################################
  
  ## If IP is unreplicated:
  if [[ $(ls ${BASENAME_IP}*dedup.bam | wc -l) < 2 ]]; then
    echo '[ERROR]: PePr requires at least n=2 for ChIP samples.
    exit 1
    fi
  
  ## if replicates for IP:  
  CHIP=$(ls ${BASENAME_IP}*dedup.bam | xargs | awk '{gsub(" ", ",");print}')
  
  ######################################################################################################
    
  ## if replicates for IgG:
  if [[ $(ls ${BASENAME_IGG}*dedup.bam | wc -l) > 1 ]]; then
    IGG=$(ls ${BASENAME_IGG}*dedup.bam | xargs | awk '{gsub(" ", ",");print}')
    fi
    
  ## if unreplicates IgG:    
  if [[ $(ls ${BASENAME_IGG}*dedup.bam | wc -l) == 1 ]]; then
    IGG=$(ls ${BASENAME_IGG}*dedup.bam)  
    fi
  
  ## if no IgG:
  if [[ $(ls ${BASENAME_IGG}*dedup.bam | wc -l) == 0 ]]; then
    echo '[ERROR]: No IgG control found'
    exit 1
    fi
    
  ######################################################################################################  
  
  PePr \
    -c ${CHIP} \
    -i ${IGG} \
    -f ${BAMTYPE} \
    -n ${BASENAME_IP} \
    --threshold=${THRESH} \
    --peaktype=${PEAKTYPE} \
    --num-processors=2 \
    --output-directory=./PepperDir
  
  PePr-postprocess \
    --peak=./PepperDir/${BASENAME_IP}__PePr_peaks.bed \
    --chip=${CHIP} \
    --input=${IGG} \
    --file-type=${BAMTYPE} \
    --remove-artifacts  
    
  mv ./PepperDir/${BASENAME_IP}__PePr_peaks.bed.passed ./PepperDir/${BASENAME_IP}_PePr_peaks_passed.bed

}; export -f PEPPER 

## Template: PEPPER chipbase iggname bamtype thresh peaktype 

## grouplevel undiff:
#ls *_dedup.bam | grep -v 'IgG' | grep '_undiff' | awk -F "_rep" '{print $1 | "sort -u" }' | \
#  parallel "PEPPER {} IgG_undiff bam 0.01 sharp 2> {}_pepr.log"
  
## grouplevel diff:
#ls *_dedup.bam | grep -v 'IgG' | grep '_diff' | awk -F "_rep" '{print $1 | "sort -u" }' | \
#  parallel "PEPPER {} IgG_diff bam 0.01 sharp 2> {}_pepr.log"
