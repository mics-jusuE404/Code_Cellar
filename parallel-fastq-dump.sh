#!/bin/bash

#######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=71
#SBATCH --partition=normal
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=fastq-dump
#######

ls *.sra | parallel -j 2 "parallel-fastq-dump -t 36 -O ./ --tmpdir ./ -s {} --gzip --split-3"
