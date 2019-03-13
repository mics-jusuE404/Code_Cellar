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

## Use Salmon on FASTQ files, providing the basename:
IDX="/scratch/tmp/a_toen03/Genomes/mm10/Gencode_M20/salmonIDX_Gencode_M20_k31"

###################################################################################################################

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

function SALMON {
  
  if [[ ! -e ${1}_1.fastq.gz || ! -e ${1}_2.fastq.gz ]]; then
    echo '[ERROR]: At least on of the input files is missing for' $1 && exit 1
    fi
  
  salmon quant \
    -l A -i $2 -p 8 \
    --no-version-check \
    --validateMappings \
    --maxMMPExtension 7 \
    --seqBias \
    --gcBias \
    -o ${1}_salmonK31 -1 ${1}_1.fastq.gz -2 ${1}_2.fastq.gz
  
}; export -f SALMON

###################################################################################################################

## Run salmon:
if [[ "$(ls *.fastq.gz 2>/dev/null | wc -l)" == 0 ]]; then
  echo '[ERROR]: No fastq.gz files present in $pwd -- exiting' && exit 1
  fi
  
ls *_1.fastq.gz | awk -F "_1" '{print $1}' | parallel -j 8 "SALMON {} $IDX"
