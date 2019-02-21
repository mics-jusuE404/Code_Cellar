#!/bin/bash

## Quantification of chromium scRNA-seq data with salmon/alevin to Gencode reference transcriptome:

###################################################################################################################

TX2GENE="/scratch/tmp/a_toen03/Genomes/mm10/Gencode_M20/gencode.vM20.Tx2Gene.txt"
IDX="/scratch/tmp/a_toen03/Genomes/mm10/Gencode_M20/salmonIDX_Gencode_M20_k31"

###################################################################################################################

## Check if all required tools are in PATH:

function PathCheck {
  
  if [[ $(command -v $1 | wc -l) == 0 ]]; then 
    echo ${1} >> missing_tools.txt
    fi
  
}; export -f PathCheck

TOOLS=(salmon)

for i in $(echo ${TOOLS[*]}); do
  PathCheck $i; done
  
if [[ -e missing_tools.txt ]] && [[ $(cat missing_tools.txt | wc -l | xargs) > 0 ]]; then
  echo '[ERROR] Tools missing in PATH -- see missing_tools.txt' && exit 1
  fi

###################################################################################################################

function ALEVIN {
  
  (>&2 paste -d " " <(echo '[INFO]' 'Alevin for' $1 'started on') <(date))
  
  if [[ ! -e ${BASENAME}_1.fastq.gz ]] || [[ ! -e ${BASENAME}_2.fastq.gz ]]; then
  echo '[ERROR] At least one input files is missing -- exiting' && exit 1
    fi
    
  BASENAME=$1
  
  salmon alevin -l ISR -1 ${BASENAME}_1.fastq.gz -2 ${BASENAME}_2.fastq.gz --chromium \
    -i $2 --tgMap $3 -p 10 -o ${BASENAME}_alevin
    
  (>&2 paste -d " " <(echo '[INFO]' 'Alevin for' $1 'ended on') <(date))  
  
}; export -f ALEVIN

ls *:.... | parallel -j 7 "ALEVIN {BASENAME} $IDX $TX2GENE 2> {}.log"
  
