#!/bin/bash

#######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --partition=normal
#SBATCH --mem=88G
#SBATCH --time=160:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=HiC_lowlevel
#######

## Alignment of Hi-C data, each mate treated as single-end to avoid insert size issues.
## Runs two of the below functions in parallel per himem-normal node:

######################################################################################################################################

## Exit function if BAM file looks corrupted or is missing:
function ExitBam {

  (>&2 echo '[ERROR]' $1 'looks suspicious or is empty -- exiting') && exit 1
  
}; export -f ExitBam

#####################################################################################################################################

## Check if BAM is corrupted or empty:
function BamCheck {
  
  BASENAME=${1%_*}
  samtools quickcheck -q $1 && echo '' >/dev/null || ExitBam $1
  
  ## Also check if file is not empty:
  if [[ $(samtools view $1 | head -n 1 | wc -l) < 1 ]]; then
    ExitBam $BASENAME
    fi
  
}; export -f BamCheck 

######################################################################################################################################

function Fq2Bam {

  BASENAME=$1
  
  BWA_IDX=/scratch/tmp/a_toen03/Genomes/hg38/bwa_index_noALT_withDecoy/hg38_noALT_withDecoy.fa
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'started on') <(date))
  
  if [[ ! -e ${BASENAME}.fastq.gz ]] ; then
    echo '[ERROR] Input file is missing -- exiting' && exit 1
    fi

  bwa mem -v 2 -R '@RG\tID:'${BASENAME}'_ID\tSM:'${BASENAME}'_SM\tPL:Illumina' \
    -t 34 ${BWA_IDX} -A1 -B4 -E50 -L0 ${BASENAME}.fastq.gz | \
    tee >(samtools flagstat -@ 2 - > ${BASENAME}_raw.flagstat) | \
  sambamba view -l 5 -t 8 -f bam -S -o ${BASENAME}_raw.bam
    
  BamCheck ${BASENAME}_sorted.bam
    
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'ended on') <(date))
    
  }; export -f Fq2Bam

######################################################################################################################################

ls *.fastq.gz | awk -F ".fastq" '{print $1}' | parallel -j 2 "Fq2Bam {} 2>> {}.log"    
