#!/bin/bash

## Subsample BAM to a given number of reads:

function SubSample {

## see also: http://crazyhottommy.blogspot.com/2016/05/downsampling-for-bam-files-to-certain.html
FACTOR=$(samtools idxstats $1 | cut -f3 | awk -v COUNT=$2 'BEGIN {total=0} {total += $1} END {print COUNT/total}')

if [[ $FACTOR > 1 ]]
  then 
  echo '[ERROR]: Requested number of reads exceeds total read count in' $1 '-- exiting' && exit 1
fi

sambamba view -s $FACTOR -t 2 -f bam -l 5 $1

}

export -f SubSample

ls *.bam | parallel "SubSample {} 50000000 > {.}_subsampled.bam"
