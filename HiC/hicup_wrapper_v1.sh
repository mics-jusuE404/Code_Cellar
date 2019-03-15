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
THREADS=32

#####################################################################################################################################

if [[ ${GENOME} == "hg19" ]]; then
  IDX="/scratch/tmp/a_toen03/Genomes/hg19/bowtie2_idx/hg19"
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
  THREADS=$4  
  DIGEST=$5
  
  (>&2 paste -d " " <(echo '[INFO] HiCUP pipeline for' "${BASENAME}" 'started on') <(date))
  
  if [[ ! -e ${BASENAME}_HiCUP ]]; then mkdir ${BASENAME}_HiCUP; fi
  
  hicup \
    --bowtie2 ${BOWTIE2} \
    --digest ${DIGEST} \
    --index ${IDX} \
    --outdir ./${BASENAME}_HiCUP \
    --temp ./${BASENAME}_HiCUP \
    --threads ${THREADS} \
    --zip ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz 
  
  (>&2 paste -d " " <(echo '[INFO] HiCUP pipeline for' "${BASENAME}" 'ended on') <(date))
  
}; export -f HICUP

#####################################################################################################################################

ls *_1.fastq.gz | awk -F "_1.fastq.gz" '{print $1}' | \
   parallel -j 2 "HICUP {} ${BOWTIE2} ${IDX} ${THREADS} ${DIGESTFILE} 2> {}.log"
