#!/bin/bash

#### Take a BAM file, e.g. from ATAC-seq and remove reads mapped to chrM, sending the BAM to stdout
#### Assumes samtools to be in PATH:

if [[ ! -e ${1}.bai ]];
  then
  echo '[INFO]:' $i 'is not indexed. Creating it now'
  sambamba index -t 8 $1
  fi

echo '[INFO]: Now removing chrM reads'
samtools idxstats $1 | \
  cut -f 1 | \
  grep -v 'chrM' | \
  xargs samtools view -@ 8 -b $1 | \

echo '[INFO]: Finished'
