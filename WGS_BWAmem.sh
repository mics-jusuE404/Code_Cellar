#!/bin/bash

## WGS alignment pipeline 04/18 using BWA mem:
## Script tested on Intel Xeon Nodes with 32 and 64 cores and > 100GB RAM.

########################################################################################################################
########################################################################################################################

## Paths:
BWA=/home/a/a_toen03/software/bwa-0.7.16a/bwa
BWA_HG38=/scratch/tmp/a_toen03/Genomes/hg38/bwa_index_noALT_withDecoy/hg38_noALT_withDecoy.fa
SAMTOOLS=$HOME/anaconda2/bin/samtools
SAMBAMBA=/home/a/a_toen03/software/sambamba
SKEWER=$HOME/software/skewer/skewer
SAMBLASTER=$HOME/software/samblaster/samblaster
PARALLEL=$HOME/software/parallel

########################################################################################################################
########################################################################################################################

## Control if paired fastq is present:
if [[ ! -e ${BASENAME}_1.fastq.gz ]] || [[ ! -e ${BASENAME}_2.fastq.gz ]]
  then
  echo '[ERROR]: Specified files are either not present, or not suffixed with *.fastq.gz'
  exit 1
  fi

## Control if output file is already present:
if [[ -e ${BASENAME}_SortedRmdup.bam ]]; then echo $BASENAME 'already exists -- exiting' && exit 1; fi

########################################################################################################################
########################################################################################################################

echo '###################################################################################################################'
echo '[INFO]: Started on date:' && date && echo ''
echo '[MAIN]: IN PROGRESS:' $BASENAME && echo 

## Adapter- and quality trimming with skewer: standard Illumina adapter, average and trailing base quality of 25:
if [[ ! -e ${BASENAME}-trimmed-pair1.fastq ]]
  then
  echo '[MAIN]: Adapter/Quality trimming'
  $SKEWER --quiet -n -q 25 -Q 25 -o $BASENAME -t 8 -m pe -l 25 ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz
  echo ''
  fi

## Align with BWA mem, mark dups with SAMblaster:
echo '[MAIN]: Alignment, sorting, markdup and samblaster (disco, splitters)'
$BWA mem -R '@RG\tID:'${BASENAME}'_ID\tSM:'${BASENAME}'_SM' -v 2 -t 24 ${BWA_HG38} ${BASENAME}-pair1.fastq ${BASENAME}-pair2.fastq | \
  $SAMTOOLS fixmate -@ 4 -O SAM - - | \
  $SAMBLASTER --addMateTags --ignoreUnmated -e -d ${BASENAME}_discordant.sam -s ${BASENAME}_splitter.sam -u ${BASENAME}_unmapped.fastq | \
  $SAMBAMBA view -f bam -S -l 0 -t 4 -o /dev/stdout /dev/stdin | \
  $SAMBAMBA sort -m 30G --tmpdir=./ -l 6 -t 32 -o ${BASENAME}_SortedRmdup.bam /dev/stdin
 
echo '[INFO]: Finished alignment on date:'
date && rm ./${BASENAME}*trim*.fastq

echo '[MAIN]: Sorting SAM files from samblaster'
ls ${BASENAME}*.sam | \
  $PARALLEL "$SAMTOOLS sort -m 2G -@ 8 -o {.}.bam {}"
echo '###################################################################################################################'

