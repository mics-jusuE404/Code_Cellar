#!/bin/bash

## Take an indexed BAM file and create a CPM-normalized bigwig, using
## - mosdepth for a per-base bed.gz (=bedgraph)
## - mawk for custom to-CPM conversion and
## bedGraphToBigWig (kentutils) for bigwig conversion

######################################################################################################
######################################################################################################

BAM=$1
CHROMSIZES=$2

## 1. Make the depth file:
mosdepth -t 8 ${1%.bam} $1

## 2. Calculate scaling factor:
SCALE_FACTOR=$(bc <<< "scale=8;1000000/$(samtools idxstats $1 | awk '{SUM+=$3} END {print SUM}')") 

## 3. Normalize and write to compressed bedGraph (can be read by IGV):
mawk -v SF=${SCALE_FACTOR} 'OFS="\t" {print $1, $2, $3, $4*SF}' <(bgzip -c -d -@ 8 ${BAM%.bam}.per-base.bed.gz) | \
  sort -k1,1 -k2,2n --parallel=8 | bgzip -@ 4 > ${BAM%.bam}_CPM.bedGraph.gz
rm ${BAM%.bam}.per-base
rm ${BAM%.bam}.*mosdepth*
