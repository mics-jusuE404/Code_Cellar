#!/bin/bash

#### Usage: ./sort_generic in_file vcf > out_file

## Uses mawk for increased speed in comparison to standard awk:

FILE=$1
TYPE=$2

if [[ $TYPE == "vcf" ]]
  then
  cat $FILE | mawk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}'
  fi
 
if [[ $TYPE == "gtf" ]]
  then
  cat $FILE | mawk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}'
  fi
  
if [[ $TYPE != "vcf" ]] || [[ $TYPE != "gtf" ]]
  then
  echo '$2 must be vcf or gtf -- exiting'
  exit 1
  fi
