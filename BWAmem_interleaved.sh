#!/bin/bash

## Saw recently that cutadapt (unlike skewer) can accept and write interleaved fq and pipe that to BWA,
## probably saving time and disk space when files are large.
## Note: When using multithreaded cutadapt, one must use the python3 version!


## add option for optional single-end input

BASENAME=$1
ADAPTER1=$2
ADAPTER2=$3
BWA_IDX=/path/to/(...)


PAIRED_2_INTERLEAVED="seqtk mergepe ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz"
BWAMEM_INTERLEAVED="bwa mem -v 2 -p ${BWA_IDX}"  
  
${PAIRED_2_INTERLEAVED} | \
  cutadapt -a $ADAPTER1 -A $ADAPTER2 --interleaved -m 18 --max-n 0.1 - | \
  ${BWAMEM_INTERLEAVED} /dev/stdin | \
  samtools fixmate -m - - | \
  samblaster --ignoreUnmated --addMateTags | \
  sambamba view -S -f bam -t 1 -l 0 /dev/stdin | \
  sambamba sort -m 5G -l 5 -t 1 --tmpdir=./ -o ${BASENAME}_raw.bam /dev/stdin 
  
exit  
  
