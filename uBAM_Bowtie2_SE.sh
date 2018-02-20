#!/bin/bash

## Simple generic alignment and filtering pipeline for standard Illumina single-end reads, e.g. ChIP-seq,
## given that data are provided as uBAM file.
## Requires skewer, samtools, sambamba and bowtie2 in PATH
## Requires that the bowtie2 index files are provided below as BOWTIE2_IDX
## Assumes to be in the same dir as the uBAMs
## Written by Alexander Toenges (a.toenges@uni-muenster.de) 2017

#GENOMES AND CHROMSIZE FILES:
BOWTIE2_IDX="/path/to/idx"

#################################################################################################################################
#################################################################################################################################
  
## Start pipeline:
for i in *.bam 
  do 
  BASENAME=${i%.bam}
  echo '[START]:' $i && echo ''
  
  ## Write uBAM to fastq, trim for base quality and degenerated reads and write as uncompressed fastq to disk:
  echo '[MAIN]: Bam2Fq & Adapter/Quality Trimming' && echo ''
  samtools fastq -@ 16 $i | \
    skewer --quiet -Q 25 -q 30 -n -t 16 -o ${BASENAME} -
  
  ## Align with bowtie2 and sort BAM file:
  echo '[MAIN]: Alignment & Sorting'
  bowtie2 -x ${BOWTIE2_IDX} -U ${BASENAME}-trimmed.fastq --quiet --very-sensitive -p 16 | \
  sambamba view -S -f bam -l 0 -t 2 /dev/stdin | \
  sambamba sort -l 9 -t 16 --tmpdir=./ -t 6 -m 20G -o ${BASENAME}_raw.bam /dev/stdin
  sambamba flagstat -t 6 ${BASENAME}_raw.bam > ${BASENAME}_raw.flagstat
  rm ${BASENAME}-trimmed.fastq
  
  ## Filter out unmapped reads and MAPQ below 30, kick out duplicates:
  echo '[MAIN]: Quality Filter & Rmdup'
  samtools view -bhu -f 0 -q 30 -@ 2 ${BASENAME}_raw.bam | \
  samtools rmdup -s - ${BASENAME}_rmdup.bam 
  sambamba index -t 6 ${BASENAME}_rmdup.bam
  sambamba flagstat -t 6 ${BASENAME}_rmdup.bam > ${BASENAME}_rmdup.flagstat
  
  done

echo 'Finished #################################################################'  
