#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --partition=normal
#SBATCH --mem=70G
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=bam2bigwig_TMM
#SBATCH --output=bam2bigwig_TMM.log

## Submit bam2bigwig_TMM_v1.sh to SLURM:

./bam2bigwig_TMM_v1.sh \
-b "*.bam" \
--atacseq \
--cutsites \
--peaks <peaks.bed> \
--mean \
--jobs 4 \
--threads 16
