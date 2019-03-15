#!/bin/bash

## Wrapper script to run the HiCUP pipeline on all fastq file pairs in $pwd.
## Written by Alexander Toenges (a.toenges@uni-muenster.de) -- 2019

#####################################################################################################################################

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --partition=normal
#SBATCH --mem=80G
#SBATCH --time=80:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=HiCUP

#####################################################################################################################################

GENOME="hg19"
CUTTER="HindIII"
BOWTIE2="/home/a/a_toen03/anaconda3_things/anaconda3/bin/bowtie2"

#####################################################################################################################################

if [[ ${GENOME} == "hg19" ]]; then
  IDX="/scratch/tmp/a_toen03/Genomes/hg19/bowtie2_idx"
  fi
  
if [[ ${CUTTER} == "HindIII" ]]; then
  DIGESTFILE="/scratch/tmp/a_toen03/Genomes/hg19/hicup_digested/Digest_hg19_HindIII_None_13-39-56_15-03-2019.txt.gz"
  fi  
  
if [[ ! -e $DIGESTFILE ]]; then
  echo '[ERROR]: Digestion file does not exist -- exit' && exit 1
  fi
  
##################################################################################################################################### 
  
## Check if required tools are in PATH:
if [[ -e missing_tools.txt ]]; then rm missing_tools.txt; fi

function PathCheck {
  
  if [[ $(command -v $1 | wc -l) == 0 ]]; then 
    echo ${1} >> missing_tools.txt
    fi
  
}; export -f PathCheck

TOOLS=(bowtie2 samtools hicup hicup_deduplicator hicup_digester hicup_filter hicup_mapper hicup_truncater)

for i in $(echo ${TOOLS[*]}); do
  PathCheck $i; done
  
if [[ -e missing_tools.txt ]] && [[ $(cat missing_tools.txt | wc -l | xargs) > 0 ]]; then
  echo '[ERROR] Tools missing in PATH -- see missing_tools.txt' && exit 1
fi  

#####################################################################################################################################

## Wrapper function for HiCUP:
function HICUP {

  BASENAME=$1
  BOWTIE2=$2
  IDX=$3
  OUTDIR=$4
  THREADS=$5  
  DIGEST=$6
  
  (>&2 paste -d " " <(echo [INFO] HiCUP pipeline for' "${BASENAME}" 'started on') <(date))
  
  mkdir ${BASENAME}_HiCUP
  
  hicup \
    --bowtie2 ${BOWTIE2} \
    --digest ${DIGEST} \
    --format Sanger \
    --index ${IDX} \
    --outdir ${OUTDIR} \
    --temp (tmp) \
    --threads ${THREADS} \
    --zip \
    ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz 
  
  (>&2 paste -d " " <(echo [INFO] HiCUP pipeline for' "${BASENAME}" 'ended on') <(date))
  
}; export -f HICUP

#####################################################################################################################################

ls *_1.fastq.gz | awk -F "_" '{print $1}' | \
   parallel -j 2 "HICUP {} ${BOWTIE2} ${IDX} {}_HiCUP ${THREADS} ${DIGESTFILE} 2> {}.log"
