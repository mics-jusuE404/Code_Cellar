#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --partition=normal
#SBATCH --mem=100G
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=DNAseq_lowlevel_v3
#SBATCH --output=DNAseq_lowlevel_v3.log

## Submit DNAseq_lowlevel_v3.sh to SLURM

./DNAseq_lowlevel_v3.sh \
--genome mm10 \
--atacseq \
--jobs 4 \
--threads 16 \
--seqtype fq_paired \
--macsPeaks

