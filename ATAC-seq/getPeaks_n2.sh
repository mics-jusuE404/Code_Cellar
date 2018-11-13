#!/bin/bash

## Assumung a duplicate of ATAC-seq experiments, call peaks on all datasets present in $(pwd),
## assuming that for a given cell type, the full BAM and one with only the read pairs with isize < 100bp are present.
## Then only keeps peaks that are present in both replicates of $full, removes blacklisted regions (output: *full.bed)
## and then also intersects the reproducible 100bp peaks with this *full.bed, keeping only those 100bp regions,
## that overlap with full, output *_NFR.bed
## The *NFR.bed is a conservative list of high-confidence NFR regions.

BLACKLIST_HG38="/scratch/tmp/a_toen03/Genomes/hg38/Blacklists/hg38_consensusBL.bed"

source activate py27

## Call peaks on full and 100bp dataset:
ls *_sorted.bam | awk -F "_sorted.bam" '{print $1}' | \
  parallel "macs2 callpeak -t {}_sorted.bam -n {} -g hs --nomodel -f BAMPE -q 0.01 --keep-dup=all --tempdir ./"

## Intersect full dataset:
ls *rep*_peaks.narrowPeak | awk -F "_rep" '{print $1}' | sort -k1,1 -u | \
  parallel "bedtools intersect -a {}_rep1_peaks.narrowPeak -b {}_rep2_peaks.narrowPeak | \
  bedtools intersect -v -a - -b /scratch/tmp/a_toen03/Genomes/hg38/Blacklists/hg38_consensusBL.bed | cut -f1-3 > {}_full.bed"

## intersect the 100bp replicates, using the overlapping coordinates and keep if they overlap with a peak of the full dataset:
ls *_rep*_sorted_isize100_peaks.narrowPeak | awk -F "_rep" '{print $1}' | sort -k1,1 -u | \
  parallel "bedtools intersect -a {}_rep1_sorted_isize100_peaks.narrowPeak -b {}_rep2_sorted_isize100_peaks.narrowPeak | \
  cut -f1-3 | bedtools intersect -wa -a - -b {}_full.bed > {}_NFR.bed"
