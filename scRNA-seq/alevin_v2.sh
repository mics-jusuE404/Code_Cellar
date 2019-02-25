#!/bin/bash

#######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --partition=normal
#SBATCH --mem=80G
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=salmon_TXquant
#SBATCH --output=salmon_quant.log
#######

###################################################################################################################

##
## Quantification of chromium scRNA-seq data with salmon/alevin to Gencode reference transcriptome.
## Alevin is a (quote from abstract of preprint):
## => a fast end-to-end pipeline to process droplet-based single cell RNA sequencing data, 
##    which performs cell barcode detection, read mapping, unique molecular identifier deduplication, 
##    gene count estimation, and cell barcode whitelisting
##

###################################################################################################################

TX2GENE="/scratch/tmp/a_toen03/Genomes/mm10/Gencode_M20/scRNA_stuff/gencode.vM20.Tx2Gene.txt"
MITORNA="/scratch/tmp/a_toen03/Genomes/mm10/Gencode_M20/scRNA_stuff/gencode.vM20_mtGenes.txt"
RRNA="/scratch/tmp/a_toen03/Genomes/mm10/Gencode_M20/scRNA_stuff/gencode.vM20_rrnaGenes.txt"
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

## Main Alevin function for 
function ALEVIN {
  
  BASENAME=$1
  (>&2 paste -d " " <(echo '[INFO]' 'Alevin for' $1 'started on') <(date))
  
  if [[ ! -e ${BASENAME}_1.fastq.gz ]] || [[ ! -e ${BASENAME}_2.fastq.gz ]]; then
  echo '[ERROR] At least one input files is missing -- exiting' && exit 1
    fi
    
  BASENAME=$1
  
  salmon alevin -l ISR -1 ${BASENAME}_1.fastq.gz -2 ${BASENAME}_2.fastq.gz --chromium \
    -i $2 --tgMap $3 --mrna $4 --rrna $5 -p 10 -o ${BASENAME}_alevin
    
  (>&2 paste -d " " <(echo '[INFO]' 'Alevin for' $1 'ended on') <(date))  
  
}; export -f ALEVIN

###################################################################################################################

ls *_1.fastq.gz | awk -F "_1.fastq.gz" '{print $1}' | \
  parallel -j 7 "ALEVIN {} $IDX $TX2GENE $MITORNA $RRNA 2> {}.log"
