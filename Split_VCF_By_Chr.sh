#!/bin/bash

## Assumes uncompressed VCF, extracts entries by chromosomes using tabix and parallel:
if [[ ${1: -4} != ".vcf" ]]
  then
  echo '[ERROR]: Please provide an uncompressed vcf file'
  exit; fi

## Get chromosomes that are present in the VCF:
cat $1 | mawk '$1 ~ /^#/ {next} {print $1 | "sort -k1,1 -u"}' > ${1%.vcf}_chrs.txt

## Compress & index file with tabix:
cat $1 | mawk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' | bgzip -c > ${1%.vcf}_sorted.vcf.gz
tabix -p vcf ${1%.vcf}_sorted.vcf.gz

## Extract files by chr:
cat ${1%.vcf}_chrs.txt | parallel "tabix -h ${1%.vcf}_sorted.vcf.gz {} > ${1%.vcf}_{}.vcf"

## Usage: ./split_vcf_by_chr.sh input.vcf
