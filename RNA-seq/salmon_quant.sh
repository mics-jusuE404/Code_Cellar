#!/bin/bash
#######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=70
#SBATCH --partition=hims
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --output=salmon_quant.log
#######

## Use Salmon on paired FASTQ files, providing the basename:

IDX=/scratch/tmp/a_toen03/Genomes/hg38/Gencode_v28/gencode.v28.salmonIDX31

## libtype ISR so first strand, stranded, paired-end
ls *_1.fastq.gz | awk -F "_1" '{print $1}' | \
  parallel "salmon quant -l A -i $IDX -p 8 --validateMappings --seqBias --gcBias -o {}_salmonQuant_k31 -1 {}_1.fastq.gz -2 {}_2.fastq.gz"
