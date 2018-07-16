#!/bin/bash

## Create a count matrix with featureCounts from a single-end BAM file,
## elongating the reads to average fragment size using mawk and bedtools:
## $1 = in.bam, $2 = chromSizes.txt, $3 = in.saf, $4 = out.matrix

bedtools bamtobed -i $1 | \
  mawk 'OFS="\t" {if ($6 == "+") {$3 = $2 + 500; print $0} else if ($6 == "-") {$2 = $3 - 500; print}}' | \
  bedtools bedtobam -i - -g $2 | \
    featureCounts -a $3 -F SAF -o $4 -T 8
    
exit
