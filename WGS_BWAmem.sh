#!/bin/bash

## WGS alignment pipeline 03/18 using BWA mem:
## Script tested on Intel Xeon Nodes with 32 and 64 cores and > 100GB RAM.

######################################################################################################
######################################################################################################

# set the number of nodes
#SBATCH --nodes=1

# set the number of CPU cores per node
#SBATCH --ntasks-per-node=64

# set a partition
#SBATCH --partition=normal

# set max wallclock time
#SBATCH --time=48:00:00 

#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=a.toenges@uni-muenster.de

## Basename must be added using an external command concatenating this upper constant part of the script,
## then the basename, then the actual pipeline

BASENAME=SRRXXXXXX

######################################################################################################
######################################################################################################

## Paths:
BWA=/home/a/a_toen03/u0dawin_compiled/bwa-0.7.16a/bwa
BWA_HG38=/scratch/tmp/a_toen03/Genomes/hg38/bwa_index_noALT_withDecoy/hg38_noALT_withDecoy.fa
SAMTOOLS=$HOME/anaconda2/bin/samtools
SAMBAMBA=$HOME/software/sambamba
SKEWER=$HOME/software/skewer/skewer
SAMBLASTER=$HOME/software/samblaster/samblaster
PARALLEL="$HOME/bin/parallel"

######################################################################################################
######################################################################################################

# Control if fastq file is gzipped:
if [[ ! -e ${BASENAME}_1.fastq.gz ]] || [[ ! -e ${BASENAME}_2.fastq.gz ]]
  then
  echo '[ERROR]: Specified files are either not present, not GZIPed or not suffixed with *.fastq.gz'
  exit 1
  fi

echo '#################################################################################'
echo '[INFO]: Started on date:'
date && echo ''
echo '[MAIN]: IN PROGRESS:' $BASENAME && echo 

# Adapter- and quality trimming with skewer:
if [[ ! -e ${BASENAME}-trimmed-pair1.fastq ]]
  then
  echo '[MAIN]: Adapter/Quality trimming'
  $SKEWER --quiet -n -q 25 -Q 25 -o $BASENAME -t 8 -m pe -l 25 ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz
  echo ''
  else
    echo '[MAIN]: Adapter/Quality skipped'
  fi

## Align trimmed data, process with fixmate, mark duplicates and extract splitter/discos with SAMBLASTER and sort with SamBamba:

echo '[MAIN]: Alignment, sorting, markdup and samblaster (disco, splitters)'
$BWA mem -R '@RG\tID:'${BASENAME}'_ID\tSM:'${BASENAME}'_SM' -v 2 -t 24 ${BWA_HG38} ${BASENAME}-trimmed-pair1.fastq ${BASENAME}-trimmed-pair2.fastq | \
  $SAMTOOLS fixmate -@ 4 -O SAM - - | \
  $SAMBLASTER --addMateTags --ignoreUnmated -e -d ${BASENAME}_discordant.sam -s ${BASENAME}_splitters.sam -u ${BASENAME}_unmapped.fastq | \
  $SAMBAMBA view -f bam -S -l 0 -t 4 -o /dev/stdout /dev/stdin | \
  $SAMBAMBA sort -m 30G --tmpdir=./ -l 6 -t 16 -o ${BASENAME}_SortedRmdup.bam /dev/stdin
  
######################################################################################################
######################################################################################################

echo '[INFO]: Finished alignment on date:'
date
rm ./${BASENAME}*trim*.fastq
if [[ -d ./${BASENAME}_TMPSORTDIR ]]
  then
  rm -r ./${BASENAME}_TMPSORTDIR
  fi
  
######################################################################################################
######################################################################################################

## Move SAM files from Samblaster in separate folder:
if [[ ! -d SAMBLASTER_files ]]
  then
  mkdir SAMBLASTER_files
  fi
mv ${BASENAME}*.sam ./SAMBLASTER_files
mv ${BASENAME}_unmapped.fastq ./SAMBLASTER_files
