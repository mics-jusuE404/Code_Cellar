#!/bin/bash

## Make reference regulome from the Genrich calls and then a count matrix, 
## assumes to be in the ATAC-seq dir, ./GenrichDir:

## merge:
cat \
  *genrich_FDR1perc.narrowPeak
  | sort -k1,1 -k2,2n --parallel=8 -S10G \
  | bedtools merge -i - \
  | tee reference_regulome.bed \
  | awk 'OFS="\t" {print $1":"$2+1"-"$3, $1, $2+1, $3, "+"}' > reference_regulome.saf 
  
## counts:
featureCounts -a reference_regulome.saf  -F SAF --read2pos 5 -T 16 -o reference_regulome.featurecounts *_dedup.bam

## also 10kb windows across the genome for potential normalization in csaw:
bedtools makewindows -g ../tmp_chromSizes.txt -w 10000 \
  | mawk 'OFS="\t" {print $1":"$2+1"-"$3, $1, $2+1, $3, "+"}'> mm10_windows10kb.saf

featureCounts -a mm10_windows10kb.saf  -F SAF --read2pos 5 -T 16 -o mm10_windows10kb.featurecounts *_dedup.bam

for (i in *.featurecounts); do
  awk 'OFS="\t" {if(NR>1) print | "cut -f2,3,4,7-"}' $i > ${%i.featurecounts}_countmatrix.tsv
  done
