#!/bin/bash

## Use stringtie on an RNA-seq cohort to:
## 1) assembly transcripts (together with reference GTF) from each sample
## 2) merge transcripts into a consensus GTF
## 3) estimate abundance based on the consensus GTF for each sample, creating a count matrix for further analysis.

## Version v_1.0, May 2018

###################################################################################################################
###################################################################################################################

REF_GTF=$1
CPU=$2

## Assumes that script is in the same directory as the coordinate-sorted RNA-seq BAM files,
## assumes that lane replicates have already been merged into single BAM files,
## assumes that BAMs are *_sorted.bam
## assumes stringtie and GNU parallel in PATH:

if [[ $# -eq 0 ]] ; then
    echo ''
    exit 0
fi

## 1) Assemble transcripts for every sample:
echo '[MAIN]: Assembling transcripts'
find . -name \*_sorted.bam | sed 's|^./||' | \
  awk -F "_sorted.bam" '{print $1}' | \
  parallel -j $CPU "stringtie {}_sorted.bam -G $REF_GTF -p 1 -l {} -o ./{}_stringtie.gtf -A {}_gene_abund.gtf -B"
  
## Merge:
echo '[MAIN]: Merging GTF files:'
stringtie --merge -G $REF_GTF -o cohort_merged.gtf *_stringtie.gtf

## Estimate abundances based on merged GTF:
echo '[MAIN]: Estimating abundances based on merged GTF file:'
find . -name \*_sorted.bam | sed 's|^./||' | \
  awk -F "_sorted.bam" '{print $1}' | \
  parallel -j $CPU "stringtie {}_sorted.bam -G cohort_merged.gtf -p 1 -o ./{}_merged_abund.gtf"
  
## Collapse *_merged_abund.gtf into single file:  

