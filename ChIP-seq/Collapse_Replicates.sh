#!/bin/bash

## Combine all files belonging to one biological replicate into one fastq.gz file,
## given they have read length >= 36bp.

CONDITIONAL_CAT () {
  BASENAME=$1
  
  echo '[INFO]: Working on' $BASENAME
  
  while read p; do
  
    ## If the read length is shorter 36bp, baclklist file and do not cat:
    if [[ $(zcat $p | head -n 40000 | mawk '{if(NR%4==2) print length($1)}' | sort -k1,1n -u) < 36 ]]
      then
      echo ${p} >> ${BASENAME}_excludeCatList.txt
      else
        
        ## If there is only one file and not multiple, simply rename:
        if [[ $(ls ${BASENAME}*.fastq.gz | wc -l) == 1 ]]; then
          mv $(ls ${BASENAME}*.fastq.gz) ${BASENAME}_combined.fastq.gz
          
          ## As only one file, go to next iteration:
          continue
          fi
        
        ## If multiple files per biol-rep. start combining:
        cat ${p} >> ${BASENAME}_combined.fastq.gz
        echo ${p} >> ${BASENAME}_includeCatList.txt
        
      fi
      
    done < <(ls ${BASENAME}*.fastq.gz)   
}; export -f CONDITIONAL_CAT    

## Run on all files in $(pwd)
ls *fastq.gz | \
  awk -F "_Tech" '{print $1}' | sort -u | \
    parallel -j 4 "CONDITIONAL_CAT {}"

