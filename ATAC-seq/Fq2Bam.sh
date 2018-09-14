#!/bin/bash


function Fq2Bam {

  BASENAME=$1
  
  ## Nextera adapter:
  ADAPTER1="CTGTCTCTTATACACATCT"
  ADAPTER2="CTGTCTCTTATACACATCT"
  
  BWA_IDX=/scratch/tmp/a_toen03/Genomes/hg38/bwa_index_noALT_withDecoy/hg38_noALT_withDecoy.fa
  
  ####################################################################################################################################
  
  echo '#############################################################################################################################'
  echo '[START]' $BASENAME 'on:' && date && echo ''
  
  seqtk mergepe ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz | \
    cutadapt -j 4 -a $ADAPTER1 -A $ADAPTER2 --interleaved -m 18 --max-n 0.1 - | \
    bwa mem -v 2 -R '@RG\tID:'${BASENAME}'_ID\tSM:'${BASENAME}'_SM\tPL:Illumina' -p -t 16 ${BWA_IDX} /dev/stdin | \
    samtools fixmate -m -@ 2 -O SAM - - | \
    samblaster --ignoreUnmated | \
    sambamba view -f bam -S -l 1 -t 4 -o /dev/stdout /dev/stdin | \
    sambamba sort -m 200G --tmpdir=./ -l 6 -t 16 -o ${BASENAME}_raw.bam /dev/stdin
    
    samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -v 'chrM|_random|chrU|chrEBV|\*' | \
      xargs sambamba view -f bam -t 8 --num-filter=1/1028 --filter='mapping_quality > 0' \
      -o ${BASENAME}_sorted.bam ${BASENAME}_raw.bam
    
  echo '[END]' $BASENAME 'on:' && date
  echo '#############################################################################################################################'
}

export -f Fq2Bam

ls *_1.fastq.gz | awk -F "_1.fastq.gz" '{print $1}' | \
  parallel -j 8 "Fq2Bam {}"
