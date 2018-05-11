#!/bin/bash

## Simple RNA-seq pipeline, accepting paired-end fastq files, producing sorted HISAT2 alignment.
## Written by Alexander Toenges (a.toenges@uni-muenster.de) 2017
## USAGE: RNAseq_hisat2.sh <Basename>

################################################################################################################################################
################################################################################################################################################

echo '[INFO]: Started on date:' && date && echo ''

BASENAME=$1

## Executables:
SKEWER=
HISAT2=
SAMBAMBA=
FEATURECOUNTS=

## HISAT2 index, splice-site file and genome GTF:
HISAT2_IDX=
SPLICE_FILE=
GTF=

## Check if fastq files are present
if [[ ! -e ${BASENAME}_1.fastq.gz ]] || [[ ! -e ${BASENAME}_1.fastq.gz ]]; then echo '[ERROR]: Input files not present -- exiting'; fi

################################################################################################################################################
################################################################################################################################################

##: Trim adapters with skewer
# --- Illumina's old paired-end adapters ::: -x fwd -y rev
# --- remove degenerated (NNNNN) reads ::: -n
# --- trim 3' until hitting a base with qual >= 20 ::: -q 20
# --- remove reads with mean quality below 10 ::: -Q 10
# --- multithreaded ::: -t INT

## Prepare output directory:
if [[ ! -d fastq_trimmed ]]; then mkdir fastq_trimmed; fi

echo '[MAIN]: Adapter/Quality trim for sample' $1
$SKEWER -m pe -n --quiet -q 30 -Q 30 -t 8 \
   -o ./fastq_trimmed/${1} ${1}_1.fastq.gz ${1}_2.fastq.gz

################################################################################################################################################
################################################################################################################################################

cd ./fastq_trimmed   

## Check for output directory:
if [[ ! -d bam_sorted ]]; then mkdir bam_sorted; fi

## Align with HISAT2 followed by sorting and indexing with Sambamba:
echo '[MAIN]: HISAT2 for sample' $1
$HISAT2 -p 32 -X 1000 --known-splicesite-infile $SPLICE_FILE --summary-file ${1}_hisat2_report.log -x $HISAT2_IDX -1 ${1}-trimmed-pair1.fastq -2 ${1}-trimmed-pair2.fastq | \
  $SAMBAMBA view -S -f bam -l 0 -p -t 4 /dev/stdin | \
  $SAMBAMBA sort -m 20G -l 5 -t 32 -o ./bam_sorted/${1}_sorted.bam /dev/stdin && rm ${1}-trimmed*.fastq
  
################################################################################################################################################
################################################################################################################################################

cd ./bam_sorted

echo '[MAIN]: featureCounts for sample' $1 'on' $GTF
## Assign reads to GTF exons:
# -a = the GTF/GFF file
# -F = specify the format of -a 
# -p = data are paired-end
# -T = set number of threads
# -P  = only consider pairs with ISIZE defined by -d & -D, default 50-600bp
# -o  = output file
$FEATURECOUNTS -a $GTF -F GTF -p -T 8 -P -o ${BASENAME}_countMatrix.txt ${1}_sorted.bam

################################################################################################################################################
################################################################################################################################################

echo '[INFO]: Ended on date:' && date 
