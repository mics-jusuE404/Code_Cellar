#!/bin/bash

## For each ATAC-seq replicate (n=2 per cell type), call peaks individually at an FDR of 0.001,
## then merge peak lists of the pair and filter against the hg38 blacklist.
## Assumes duplicate- and MAPQ filtered BAM files of the replicate in same directory.

BLACKLIST_HG38="/scratch/tmp/a_toen03/Genomes/hg38/Blacklists/hg38_consensusBL.bed"

## Call peaks:
ls *_sorted.bam | awk -F "_sorted.bam" '{print $1}' | \
  parallel "macs2 callpeak -t {}_sorted.bam -n {} -g hs --nomodel -f BAMPE -q 0.01 --keep-dup=all"

## Merge peak lists and filter against blacklist:
ls *peaks.narrowPeak | awk -F "_rep" '{print $1}' | sort -u | \
  parallel "cat {}_rep*.narrowPeak | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - | \
    bedtools intersect -v -a - -b ${BLACKLIST_HG38} | \
    cut -f1-3 > {}_peaks_merged.bed"
