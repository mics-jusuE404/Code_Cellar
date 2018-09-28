#!/bin/bash

## Given a metadata.tsv file from encode in $(pwd), download the fastq files and rename with experimental target,
## the biological replicate number and accession number:

## Download paths from metadata:
grep 'fastq.gz' metadata.tsv | cut -f1 | while read p; do echo 'https://www.encodeproject.org/files/'${p}'/@@download/'${p}'.fastq.gz' >> download.txt;  done

## Download:
cat download.txt | head -n 2 | parallel -j 6 "wget {}"

## Rename fastq in the form EXPERIMENTTARGET_ACCESSION_BIOLREP.fastqgz, skipping unlisted files:
for i in *.fastq.gz; do 
  if [[ $(grep $i metadata.tsv) == "" ]]; then continue; fi
  mv $i "$(grep $i metadata.tsv | cut -f13)_${i%.fastq.gz}_biolRep$(grep $i metadata.tsv | cut -f25).fastq.gz"
  done
