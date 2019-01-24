#!/bin/bash

## Alignment of Hi-C data, each mate treated as single-end to avoid insert size issues:

function Fq2Bam {

  BASENAME=$1
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'started on') <(date))
  
  if [[ ! -e ${BASENAME}.fastq.gz ]] ; then
    echo '[ERROR] Input file is missing -- exiting' && exit 1
    fi

  mkdir ${BASENAME}_tmpDir
  
  bwa mem -v 2 -R '@RG\tID:'${BASENAME}'_ID\tSM:'${BASENAME}'_SM\tPL:Illumina' -t 16 ${BWA_IDX} ${BASENAME}.fastq.gz | \
    samblaster --ignoreUnmated | \
    sambamba view -f bam -S -l 0 -t 4 -o /dev/stdout /dev/stdin | \
      tee >(samtools flagstat -@ 2 - > ${BASENAME}_raw.flagstat) | \
    samtools sort -l 5 -T ./${BASENAME}_tmpDir/${BASENAME}_tmp -@ 16 -m 2G | \
      tee ${BASENAME}_sorted.bam | \
    samtools index -@ 4 - ${BASENAME}_sorted.bam.bai
    
    rm -r ${BASENAME}_tmpDir
    
    BamCheck ${BASENAME}_sorted.bam
    
    }; export -f Fq2Bam
    
ls *.fastq.gz | awk -F ".fastq" '{print $1}' | parallel -j 4 "Fq2Bam {}"    
