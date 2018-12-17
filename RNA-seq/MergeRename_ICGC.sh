#!/bin/bash

## See README.md for infos

###################################################################################################################

## List of unique sample IDs:
ls *R1.fastq.gz | \
  awk -F "-tumor_" '{print $2}' | \
  awk -F "_RNA-" '{print $1}' | \
  awk '$1 == "" {next} {print | "sort -k1,1 -u"}' > IDs_unique.txt

###################################################################################################################

## Check if files in $pwd:
if [[ "$(cat IDs_unique.txt | wc -l)" == 0 ]]; then 
  echo ''
  echo '[ERROR]: No files present in $pwd' && exit 1
  echo ''
  fi

###################################################################################################################

## Merge reps into *_combined_1/2.fastq.gz if multiple runs per sample are present or simply rename file
## in case only one replicate is present:
function MergeRename {
  
  ID=$1
  
  ## Check if both numbers of R1 files match R2	files for each ID:
  if [[	"$(ls *${ID}*R1.fastq.gz 2>/dev/null | wc -l)" != "$(ls *${ID}*R2.fastq.gz 2>/dev/null | wc -l)" ]]; then
    echo '[ERROR]: For' $ID 'the numbers of R1 files are != R2 -- skipping this	ID' >> report_error.log
    exit 1
    fi

  ## If	ok, start merging if > 1 file per R1/2:
  if [[ "$(ls *${ID}*R1.fastq.gz | wc -l)" > 1 ]]; then
    echo 'Merging' $ID 'R1' >> report_merged.log
    find . -maxdepth 1 -name "*${ID}*R1.fastq.gz" | xargs cat >	${ID}_RNAseq_combined_1.fastq.gz
    echo 'Merging' $ID 'R2' >> report_merged.log
    find . -maxdepth 1 -name "*${ID}*R2.fastq.gz" | xargs cat >	${ID}_RNAseq_combined_1.fastq.gz
    echo '' >> report_merged.log
    fi
  
  ## If ok but only one file per R1/2, rename:
  if [[ "$(ls *${ID}*R1.fastq.gz | wc -l)" == 1 ]]; then
    echo 'Renaming' $ID 'R1' >> report_renamed.log
    mv *${ID}*R1.fastq.gz ${ID}_RNAseq_combined_1.fastq.gz
    echo 'Renaming' $ID 'R2' >> report_renamed.log
    mv *${ID}*R2.fastq.gz ${ID}_RNAseq_combined_2.fastq.gz
    echo '' >> report_renamed.log
    fi
    
}; export -f MergeRename

###################################################################################################################

cat IDs_unique.txt | parallel -j 6 MergeRename {}
  
