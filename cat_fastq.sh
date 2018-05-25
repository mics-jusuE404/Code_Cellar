#!/bin/bash

## Concatenates fastq files (*_1.fastq.gz, *_2.fastq.gz), taking an input list with the basenames as $1.

if [[ $# -eq 0 ]] ; then
    echo 'Usage: ./cat_fastq.sh input.list'
    exit 0
fi

while read p; do
  if [[ ! -e ${p}_1.fastq.gz ]] || [[ ! -e ${p}_2.fastq.gz ]]
    then
    echo $p 'not found -- exiting' && exit 1
  fi
done < $1

if [[ -e ${1%.txt}_2.fastq.gz ]] || [[ -e ${1%.txt}_2.fastq.gz ]]; then
  echo 'Output files already present -- exiting' && exit 1
fi

echo 'All input files present -- catting:'
cat $1
while read p; do echo ${p}_1.fastq.gz; done < $1 | xargs cat >> ${1%.txt}_1.fastq.gz &
while read p; do echo ${p}_2.fastq.gz; done < $1 | xargs cat >> ${1%.txt}_2.fastq.gz
echo '' && echo 'DONE'
