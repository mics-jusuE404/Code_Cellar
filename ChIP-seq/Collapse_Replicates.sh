#!/bin/bash

## Combine all files belonging to one biological replicate into one fastq.gz file,
## given they have read length >= 36bp.

CONDITIONAL_CAT () {
  BASENAME=$1
  while read p; do
    if [[ $(zcat $p | head -n 40000 | mawk '{if(NR%4==2) print length($1)}' | sort -k1,1n -u) < 36 ]]
      then
      echo ${p} >> ${BASENAME}_excludeCatList.txt
      else
        cat ${p} >> ${BASENAME}_combined.fastq.gz
        echo ${p} >> ${BASENAME}_includeCatList.txt
      fi
    done < <(ls ${BASENAME}*.fastq.gz)   
}; export -f CONDITIONAL_CAT    

## Run on all files in $(pwd)
ls *fastq.gz | \
  awk -F "_Tech" '{print $1}' | sort -u | \
    parallel -j 4 "CONDITIONAL_CAT {}"

