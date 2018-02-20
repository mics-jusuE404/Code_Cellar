#!/bin/bash

#### Take a BAM file, e.g. from ATAC-seq and remove reads mapped to chrM
#### Assumes samtools to be in PATH:

if [[ ! -e ${1}.bai ]];
  then
  echo '[INFO]:' $i 'is not indexed. Creating it now'
  fi
samtools index $1 && echo '[INFO]:' $i 'is now indexed'

echo '[INFO]: Now removing chrM reads'
samtools idxstats $1 | \
  cut -f 1 | \
  grep -v 'chrM' | \
  xargs samtools view -hb -o ${1}_exMito.bam $1 | \

echo '[INFO]: Finished'
