#!/bin/bash

######################################################################################################################################

## ChIP-seq - alignment of single-end fastq files to hg38 assuming very short reads from ENCODE:
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
  
  ALN_IDX=/scratch/tmp/a_toen03/Genomes/hg38/bowtie_index_noALT_withDecoy/hg38_noALT_withDecoy
  
  ####################################################################################################################################
  
  echo '#############################################################' >> ${BASENAME}.log
  echo '[START]' $BASENAME 'on:' >> ${BASENAME}.log && date >> ${BASENAME}.log && echo '' >> ${BASENAME}.log
  
  ## Map with bowtie trimming to 36bp:
  seqtk trimfq -L 36 ${BASENAME}.fastq.gz | \
    bowtie --sam --quiet -m 1 --best --strata --threads 16 ${ALN_IDX2} /dev/stdin >> ${BASENAME}.log | \
    samblaster --ignoreUnmated 2>> ${BASENAME}.log | \
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
  
## Fastqc:
ls *_raw.bam | parallel "fastqc {}" 
