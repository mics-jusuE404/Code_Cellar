#!/bin/bash

######################################################################################################################################

## ChIP-seq - alignment of single-end fastq files to hg38 :::
## Assumes script in same dir as fastqs, bwa - samtools - samblaster - deeptools - fastqc
## Last update: 20.11.18

######################################################################################################################################

#######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=70
#SBATCH --partition=hims
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=ChIPseq_Align
#SBATCH --output=Fq2Bam.log
#######

function Fq2Bam {

  BASENAME=$1
 
  if [[ ! -e ${BASENAME}.fastq.gz ]]; then
    echo '[ERROR] Input file is missing -- exiting' && exit 1
    fi
  
  BWA_IDX=/scratch/tmp/a_toen03/Genomes/hg38/bwa_index_noALT_withDecoy/hg38_noALT_withDecoy.fa
  
  ####################################################################################################################################
  
  echo '#############################################################' >> ${BASENAME}.log
  echo '[START]' $BASENAME 'on:' >> ${BASENAME}.log && date >> ${BASENAME}.log && echo '' >> ${BASENAME}.log
  
  ## Map BWA mem regardless of read length:
  bwa mem -v 2 -R '@RG\tID:'${BASENAME}'_ID\tSM:'${BASENAME}'_SM\tPL:Illumina' -t 16 ${BWA_IDX} ${BASENAME}.fastq.gz 2>> ${BASENAME}.log | \
    samblaster --ignoreUnmated 2>> ${BASENAME}.log 2>> ${BASENAME}.log | \
    sambamba view -f bam -S -l 1 -t 4 -o /dev/stdout /dev/stdin 2>> ${BASENAME}.log | \
    sambamba sort -m 4G --tmpdir=./ -l 6 -t 16 -o ${BASENAME}_raw.bam /dev/stdin 2>> ${BASENAME}.log
  
  ## Remove unmapped and duplicated reads (1028):
  samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -vE 'chrM|_random|chrU|chrEBV|\*' 2>> ${BASENAME}.log | \
    xargs sambamba view -f bam -t 8 --num-filter=0/1028 --filter='mapping_quality > 19' 2>> ${BASENAME}.log \
    -o ${BASENAME}_sorted.bam ${BASENAME}_raw.bam 2>> ${BASENAME}.log 
    
  ls ${BASENAME}*.bam | parallel "sambamba flagstat -t 8 {} > {.}.flagstat 2>> ${BASENAME}.log"
    
  echo '[END]' $BASENAME 'on:' >> ${BASENAME}.log && date >> ${BASENAME}.log
  
  echo '#############################################################' >> ${BASENAME}.log
}

export -f Fq2Bam

## Alignment:
ls *.fastq.gz | awk -F ".fastq.gz" '{print $1}' | \
  parallel -j 4 "Fq2Bam {}"
  
## Bigwigs:
ls *_sorted.bam | parallel -j 4 "bamCoverage -e 400 --normalizeUsing CPM -bs 1 --bam {} -o {.}_CPM.bigwig -p 16"

## Fastqc:
ls *_raw.bam | parallel "fastqc {}" 
