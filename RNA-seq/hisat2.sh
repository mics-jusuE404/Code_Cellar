#!/bin/bash

## Simple RNA-seq pipeline, accepting paired-end fastq files, producing sorted HISAT2 alignment.
## Written by Alexander Toenges (a.toenges@uni-muenster.de) 2017
## USAGE: RNAseq_hisat2.sh <Basename>

################################################################################################################################################
################################################################################################################################################

echo '[INFO]' $1 'Started on date:' && date && echo ''

BASENAME=$1

## HISAT2 index, splice-site file and genome GTF:
HISAT2_IDX="/scratch/tmp/a_toen03/Genomes/hg38/hisat2_IDX/hg38_noALT_withDecoy"
SPLICE_FILE="/scratch/tmp/a_toen03/Genomes/hg38/hisat2_IDX/gencode.v27.annotation_SpliceSites.txt"

## Check if fastq files are present
if [[ ! -e ${BASENAME}_1.fastq.gz ]] || [[ ! -e ${BASENAME}_1.fastq.gz ]]; then echo '[ERROR]: Input files not present -- exiting'; fi

################################################################################################################################################
################################################################################################################################################

##: Trim adapters with skewer
# --- Illumina's paired-end adapters ::: -x fwd -y rev
# --- remove degenerated (NNNNN) reads ::: -n
# --- trim 3' until hitting a base with qual >= 20 ::: -q 20
# --- remove reads with mean quality below 20 ::: -Q 20
# --- multithreaded ::: -t INT

## Prepare output directory:
echo '[MAIN]: Adapter/Quality trim for sample' $1
skewer -m pe -n --quiet -q 20 -Q 20 -t 8 -o ${1} ${1}_1.fastq.gz ${1}_2.fastq.gz

################################################################################################################################################
################################################################################################################################################

## Check for output directory:
if [[ ! -d BAM ]]; then mkdir BAM; fi

## Align with HISAT2, without sorting, as mostly BAMs are only used for featureCounts quant 
## and that requires name-sorting anyway:
echo '[MAIN]: HISAT2 for sample' $1
hisat2 -p 32 -X 1000 --known-splicesite-infile $SPLICE_FILE --summary-file ${1}_hisat2_report.log -x $HISAT2_IDX -1 ${1}-trimmed-pair1.fastq -2 ${1}-trimmed-pair2.fastq | \
  samblaster | \
  sambamba view -S -f bam -l 5 -t 4 -o ./BAM/${1}_unsorted.bam /dev/stdin && rm ${1}-trimmed*.fastq
  
################################################################################################################################################
################################################################################################################################################

echo '[INFO]' $1 'Ended on date:' && date 
