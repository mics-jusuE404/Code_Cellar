#!/bin/bash

## Call peaks on ChIP-seq data downloaded from ENCODE and processed with ./Collapse_Replicates.sh.
## Keep only peaks found in both replicate

source activate py27

## => H3K4me1/2/3, H3K9ac, H3K27ac, H2A.Z and transcription factors are narrow signals
## => H3K36me3, H3K79me2, H3K27me3, H3K9me1/3 are broad signals
## One can also say that everything not broad is narrow (grep -v)

## Call:
ls *_sorted.bam | \
  grep -vE 'H3K36me3|H3K79me2|H3K27me3|H3K9me1|H3K9me3|Control' | \
  parallel "macs2 callpeak -q 0.01 --outdir ./ --tempdir ./ --verbose 0 -g hs --keep-dup=all -n {.} -t {} -c *Control*_combined_sorted.bam"

ls *_sorted.bam | \
  grep -E 'H3K36me3|H3K79me2|H3K27me3|H3K9me1|H3K9me3' | \
  parallel "macs2 callpeak --broad --broad-cutoff 0.1 --outdir ./ --tempdir ./ --verbose 0 -g hs --keep-dup=all -n {.} -t {} -c *Control*_combined_sorted.bam"

## Filter:
tba...
  
