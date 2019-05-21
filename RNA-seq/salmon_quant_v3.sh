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
MODE="PE"
TRIM="y"

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

function Trim {

  BASENAME=$1
  cutadapt -j 4 -a AGATCGGAAGAG -A AGATCGGAAGAG -m 31 --max-n 0.1 \
    -o ${BASENAME}_trimmed_1.fastq.gz \
    -p ${BASENAME}_trimmed_2.fastq.gz \
    ${BASENAME}_1.fastq.gz \
    ${BASENAME}_2.fastq.gz
    
    if [[ ! -d untrimmed ]]; then mkdir untrimmed; fi
    mv ${BASENAME}_1.fastq.gz untrimmed
    mv ${BASENAME}_2.fastq.gz untrimmed
    
}; export -f Trim

###################################################################################################################

function SALMON {
  
  if [[ ${3} == "PE" ]]; then
  
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
      -o ${1}_salmon -1 ${1}_1.fastq.gz -2 ${1}_2.fastq.gz
      
      fi
      
   if [[ ${3} == "SE" ]]; then
   
     if [[ ! -e ${1}.fastq.gz ]]; then
      echo '[ERROR]: At least on of the input files is missing for' $1 && exit 1
      fi
  
      salmon quant \
        -l A -i $2 -p 8 \
        --no-version-check \
        --validateMappings \
        --maxMMPExtension 7 \
        --seqBias \
        -o ${1}_salmon -r ${1}.fastq.gz
      
      fi
      
         
  
}; export -f SALMON

###################################################################################################################

## Optional trimming:
if [[ ${TRIM} == "y" ]]; then
  ls *_1.fastq.gz | awk -F "_1.fastq" '{print $1}' | parallel -j 12 "Trim {} 2>> trimming_report.log
  fi"
  

## Run salmon:
if [[ "$(ls *.fastq.gz 2>/dev/null | wc -l)" == 0 ]]; then
  echo '[ERROR]: No fastq.gz files present in $pwd -- exiting' && exit 1
  fi

echo "[INFO]: This is" $(salmon --version)

if [[ ${MODE} == "PE" ]]; then
  ls *_1.fastq.gz | awk -F "_1.fastq" '{print $1}' | parallel -j 8 "SALMON {} $IDX ${MODE} 2> {}_salmon.log"
  fi
  
if [[ ${MODE} == "SE" ]]; then
  ls *.fastq.gz | awk -F ".fastq" '{print $1}' | parallel -j 8 "SALMON {} $IDX ${MODE} 2> {}_salmon.log"
  fi
