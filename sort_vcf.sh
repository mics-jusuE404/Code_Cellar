#!/bin/bash

#### Usage: ./sort_vcf.sh in.vcf > out_sorted.vcf

## Uses mawk for increased speed in comparison to standard awk:
cat $1 | mawk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}'

## Can be parallelized on multiple VCF using GNU parallel, e.g:
## ls *.vcf | parallel "./sort_vcf.sh {} > {.}_sorted.vcf"
