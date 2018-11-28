#!/bin/bash
#######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=70
#SBATCH --partition=hims
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=salmon_TXquant
#SBATCH --output=salmon_quant.log
#######

## Use Salmon on FASTQ files, providing the basename:
IDX="/scratch/tmp/a_toen03/Genomes/hg38/Gencode_v28/gencode.v28.salmonIDX25"

function SALMON {
  salmon quant -l A -i $2 -p 8 --validateMappings --seqBias --gcBias -o ${1}_salmonK25 -1 ${1}_1.fastq.gz -2 ${1}_2.fastq.gz
}; export -f SALMON

## Run salmon:
ls *_1.fastq.gz | awk -F "_1." '{print $1}' | parallel -j 8 "SALMON {} $IDX"

## fastqc
ls *.fastq.gz | parallel -j 70 "fastqc {}"
