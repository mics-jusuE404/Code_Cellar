#!/bin/bash

## Write a paired-end BAM file in BED format using start of mate1 and end of mate2 as coordinates:
## Assumes samtools and bedtools in PATH:

## USAGE: ./Bam2BedPE.sh in.bam > out_bedpe.bed

if [[ $# -eq 0 ]] ; then
  echo '[ERROR]: No input file given!'
  exit 0
fi

function FUNC_QUIT {
  echo '[ERROR]:' $1 'looks corrupted or is not a BAM file -- exiting'
  exit
}

if [[ -e $1 ]]; then
  samtools quickcheck $1 && echo '[INFO]: Bam2Bed for' $1 || FUNC_QUIT
  samtools sort -l 0 -n -T ${1}_tmp.SortByNameBedpe $1 -O bam | \
    bedtools bamtobed -i - -bedpe | cut -f1,2,6,7,8,9 | sort -k1,1 -k2,2n
  else 
    echo '[ERROR]: Input file does not exist -- exiting'
    exit 1
  fi
