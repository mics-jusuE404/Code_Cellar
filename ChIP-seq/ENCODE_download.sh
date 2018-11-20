#!/bin/bash

## Script takes a metadata.tsv from ENCODE, e.g.: 
## https://www.encodeproject.org/metadata/type=Experiment&assay_term_name=ChIP-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_term_name=Karpas-422/metadata.tsv
## and downloads all fastq files via wget to $(pwd) that are either single-end or fwd-read of paired-end data.
## Then, it renames the files that are ENXXXXXX.fastq.gz into:
## Biosample_Assay_FileAccession_ExperimentTarget_Biolrep<int>_TechRep<int>.fastq.gz

########################################################################################################################################

## Main function:
LOAD_RENAME () {
  VAR_STDIN=$1
    
  if [[ $(echo "$VAR_STDIN" | cut -f29) == "paired-ended" && $(echo "$VAR_STDIN" | cut -f30) == 2 ]]; then
    exit 0
    fi
    
  ## Download file:
  wget -nv $(echo "$VAR_STDIN" | cut -f37)
    
  ## Rename:  
  mv \
    $(echo "$VAR_STDIN" | cut -f37 | awk -F "@@download/" '{print $2}') \
    $(sed 's#\-human##g' <(paste -d "" \
      <(paste -d "_" \
        <(echo "$VAR_STDIN" | cut -f7) \
        <(echo "$VAR_STDIN" | cut -f5) \
        <(echo "$VAR_STDIN" | cut -f1) \
        <(echo "$VAR_STDIN" | cut -f13) \
        <(paste -d "" <(echo 'BiolRep') <(echo "$VAR_STDIN" | cut -f25)) \
        <(paste -d "" <(echo 'TechRep') <(echo "$VAR_STDIN" | cut -f26))) \
        <(echo '.fastq.gz'))
      )
}; export -f LOAD_RENAME

########################################################################################################################################

## Run it on all metadata TSVs in $(pwd):
for i in *_metadata.tsv
  do
  echo '##################################################################'
  echo '[INFO]: Running' $i
  awk '$2 ~ /fastq/ {print}' $i | parallel -j 4 "LOAD_RENAME {}"
  echo ''
  done
  
########################################################################################################################################

