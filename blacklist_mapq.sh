#!/bin/bash

## Script takes a BAM file ($1) and outputs a BED file with all bases that have a fraction of reads with MAPQ=0 greater than
## a user-specified cutoff ($2)
## Assumes sambamba and BEDtools in PATH

## Usage: ./backlist_mapq.sh in.bam fraction_cutoff
## e.g.:  ./backlist_mapq.sh in.bam 0.2 ## to get all bases with more than 20% MAPQ=0 reads

########################################################################################################################
########################################################################################################################

if [[ $# -eq 0 ]] ; then
  echo '[ERROR]: Not enough arguments!'
  echo 'Usage: blacklist_mapq.sh in.bam mapq0_fraction[integer or float]'
  exit 0
fi

## Check if file is indexed:
if [[ ! -e ${1}.bai ]]; then
  echo '[INFO]:' $1 'not indexed -- indexing now:'
  sambamba index -t 8 $1
  fi

## Make the blacklist:
echo '[INFO]: Blacklisting all bases with MAPQ=0 fraction above' $2 '--- output file is:' ${1%.bam}_mapqBL.bed
paste \
  <(sambamba depth base -t 4 --min-coverage=0 -F 'mapping_quality == 0' --regions=chr1:1-1000 tmp.bam | cut -f1,2,3) \
  <(sambamba depth base -t 4 --min-coverage=0 -F 'mapping_quality > 0' --regions=chr1:1-1000 tmp.bam | cut -f1,2,3) | \
  awk -v FRAC="$2" 'NR>1, OFS="\t" {if (($3+1)/($6+1) >= FRAC) print $1, $2, $2+1}' | \
  bedtools merge -i - \
  > ${1%.bam}_mapqBL.bed
