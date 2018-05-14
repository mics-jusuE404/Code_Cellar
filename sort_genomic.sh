#!/bin/bash

#### Usage: ./sort_generic in_file > out_file

## Uses mawk for increased speed in comparison to standard awk:

if [[ $TYPE == "vcf" ]]
  then
  cat $1 | mawk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}'
  fi
 
if [[ $TYPE == "gtf" ]]
  then
  cat $1 | mawk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}'
  fi
  
if [[ $TYPE != "vcf" ]] || [[ $TYPE != "gtf" ]]
  then
  echo '$2 must be vcf or gtf -- exiting'
  exit 1
  fi

## Can be parallelized on multiple VCF using GNU parallel, e.g:
## ls *.vcf | parallel "./sort_vcf.sh {} > {.}_sorted.vcf"
