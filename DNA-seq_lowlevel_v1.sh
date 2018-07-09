#!/bin/bash

################################################################################################################################
################################################################################################################################

#### Script for lowlevel processing of DNA-seq data (ATAC-seq, ChIP-seq) ####
#### Outputting 1) *_raw.bam = aligned, duplicate-marked and sorted bam file including all reads.
####            2) *_sorted.bam = like 1) with multimappers and reads on mtDNA removed
####            3) *_sorted_CPM.bigwig = bigwig file from 2) normalized to CPM

################################################################################################################################
################################################################################################################################

BWA_IDX="scratch/tmp/a_toen03/Genomes/hg38/bwa_index_noALT_withDecoy/hg38_noALT_withDecoy.fa"

################################################################################################################################
################################################################################################################################

echo '[INFO]: Started on date:'
date && echo ''

## Make sure all arguments are given:
if [[ $# -ne 2 ]] ; then
  echo '[ERROR]: Too few arguments given!'
  echo 'Usage: ./DNA-seq_lowlevel_v1.sh BASENAME EXPERIMENT'
  echo 'BASENAME is the suffix of the fastq.gz file'
  echo 'EXPERIMENT can be ATAC-PE, ATAC-SE, ChIP-SE, ChIP-PE'; exit 1; fi

BASENAME=$1

if [[ $2 != "ATAC-PE" ]] && [[ $2 != "ATAC-SE" ]] && [[ $2 != "ChIP-SE" ]] && [[ $2 != "ChIP-PE" ]]; then
  echo '[ERROR]: Select the type of experiment: ATAC-PE, ATAC-SE, ChIP-SE, ChIP-PE'
  exit 1; fi

################################################################################################################################
################################################################################################################################

if [[ $2 == "ATAC-SE" ]] || [[ $2 == "ChIP-SE" ]]; then
  MODE="se"; fi
  
if [[ $2 == "ATAC-PE" ]] || [[ $2 == "ChIP-PE" ]]; then
  MODE="pe"; fi
  
## Make sure files exist:
if [[ $MODE == "pe" ]]; then
  if [[ ! -e  ${BASENAME}_1.fastq.gz ]] || [[ ! -e  ${BASENAME}_2.fastq.gz ]]; then
    echo '[ERROR]: Fastq files seem to be missing '; exit 1; fi
fi

if [[ $MODE == "se" ]]; then
  if [[ ! -e  ${BASENAME}.fastq.gz ]]; then
    echo '[ERROR]: Fastq files seem to be missing '; exit 1; fi
fi

################################################################################################################################
################################################################################################################################

SKEWER="skewer --quiet -n -q 30 -Q 25"

if [[ $2 == "ATAC-PE" ]]; then 
  SKEWER_2="-m pe -x CTGTCTCTTATACACATCT -y CTGTCTCTTATACACATCT -o ${BASENAME} ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz"; fi

if [[ $2 == "ATAC-SE" ]]; then 
  SKEWER_2="-m tail -x CTGTCTCTTATACACATCT --1 ${BASENAME}.fastq.gz"; fi
  
if [[ $2 == "ChIP-SE" ]]; then 
  SKEWER_2="-m tail -1 ${BASENAME}.fastq.gz"; fi
  
if [[ $2 == "ChIP-PE" ]]; then 
  SKEWER_2="-m pe -o ${BASENAME} ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz"; fi

################################################################################################################################
################################################################################################################################

if [[ $MODE == "se" ]]; then

  ## Single-end option:

  BWA_SE="bwa mem -v 2 -t 8 ${BWA_IDX} /dev/stdin"

  ## Trim & Align:
  $SKEWER -t 4 $SKEWER_2 | \
    ${BWA_SE}  | \
    samblaster | \
    sambamba view -S -f bam -t 1 -l 0 /dev/stdin | \
    sambamba sort -m 5G -l 5 -t 4 --tmpdir=./ -o ${BASENAME}_raw.bam /dev/stdin

  ## Filter out chrM reads, multimappers and duplicates according to FLAG from samblaster:
  samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -v 'chrM' | \
    xargs sambamba view -f bam -t 4 --num-filter=0/1028 --filter='mapping_quality > 0' \
    -o ${BASENAME}_sorted.bam ${BASENAME}_raw.bam

  fi    

################################################################################################################################
################################################################################################################################

if [[ $MODE == "pe" ]]; then

  ## Paired-end option

  BWA_PE="bwa mem -v 2 -t 8 ${BWA_IDX} ${BASENAME}-trimmed-pair1.fastq ${BASENAME}-trimmed-pair2.fastq"

  ## Trim & Align:
  $SKEWER -t 4 $SKEWER_2  
    ${BWA_PE} | \
    samblaster --ignoreUnmated | \
    sambamba view -S -f bam -t 1 -l 0 /dev/stdin | \
    sambamba sort -m 5G -l 5 -t 4 --tmpdir=./ -o ${BASENAME}_raw.bam /dev/stdin

  ## Filter out chrM reads, multimappers and duplicates according to FLAG from samblaster:
  samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -v 'chrM' | \
    xargs sambamba view -f bam -t 4 --num-filter=1/1028 --filter='mapping_quality > 0' \
    -o ${BASENAME}_sorted.bam ${BASENAME}_raw.bam
    
  fi   
  
################################################################################################################################
################################################################################################################################

## Flagstat and bigwig:
ls ${BASENAME}*.bam | parallel "sambamba flagstat -t 2 {} > {.}.flagstat"

## Browser Track (deeptools3.0 nor offers CPM normalization to avoid large numbers when using -bs 1 because 1bp = 0.001kb, which
## causes the value for each bin to be multiplied by 1/0.001 so 1000. cpm ignores the bins
BIGWIG_GENERAL="bamCoverage -p 8 --normalizeUsing CPM --bam ${BASENAME}_sorted.bam -o ${BASENAME}_sorted_CPM.bigwig -bs 1"

if [[ $MODE == "se" ]]; then
  $BIGWIG_GENERAL -e 300; fi 

if [[ $MODE == "pe" ]]; then
  $BIGWIG_GENERAL -e; fi
  
echo '[INFO]: Ended on date:'
date

################################################################################################################################
################################################################################################################################
