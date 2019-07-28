#!/bin/bash

## Call peaks with Genrich both on each sample and on sample groups indicated by BASENAME_rep*_dedup.bam

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --partition=normal
#SBATCH --mem=80G
#SBATCH --time=08:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=Genrich

## mm10 consensus blacklist for ATAC-seq:
Blacklist="/scratch/tmp/a_toen03/Genomes/mm10/mm10_consensusBL.bed"

## Sort BAM by name:
if [[ ! -e ./GenrichDir ]]; then
  mkdir GenrichDir; fi

########################################################################################################################

## Genrich requires queryname-sorted files:
ls *.bam | parallel -j 4 "samtools sort -n -@ 8 -m 2G -o ./GenrichDir/{} {}"
  
cd GenrichDir

########################################################################################################################

## Genrich for groups:
function GENRICHGROUP {
  
  BASENAME=$1
  
  ## Only if replicates are present:
  if [[ $(ls ${BASENAME}*dedup.bam | wc -l) < 2 ]]; then
    echo '[WARNING]: Only one sample found for' $BASENAME '-- skipping group-level peak calling'
    exit 0
  
  FILES=$(ls ${BASENAME}*dedup.bam | xargs | awk '{ print "\""$0"\""}')
  
  Genrich -E $2 -t $FILES -j -l 200 -q 0.01 -o ${BASENAME}_genrich_FDR1perc.narrowPeak
  
}; export -f GENRICHGROUP

ls *_dedup.bam | awk -F "_rep" '{print $1 | "sort -u"}' | parallel "GENRICHGROUP {} $Blacklist 2> {}_genrich.log"

########################################################################################################################

## Genrich for individual samples:
function GENRICHSINGLE {
  
  BASENAME=$1
  
  Genrich -E $2 -t ${BASENAME}_dedup.bam -j -l 200 -q 0.01 -o ${BASENAME}_genrich_FDR1perc.narrowPeak
  
}; export -f GENRICHSINGLE

ls *_dedup.bam | awk -F "_dedup.bam" '{print $1}' | parallel -j 8 "GENRICHSINGLE {} $Blacklist 2> {}_genrich.log"
