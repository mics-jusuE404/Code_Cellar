#!/bin/bash

## Assumes at least sorted VCF, extracts entries by chromosomes using tabix and parallel:
if [[ ${1: -3} != ".gz" ]]
  then
  echo '[INFO]: File does not seem to be compressed -- compressing now'
  bgzip $1
  tabix -p vcf ${1}.gz; fi

## Get chromosomes that are present in the VCF:
tabix ${1}.gz -l > ${1%.vcf.gz}_chr.txt && \

## Extract files by chr:
cat ${1%.vcf.gz}_chr.txt | parallel "tabix -h ${1}.gz {} > ${1%.vcf.gz}_{}.vcf"

## Usage: ./split_vcf_by_chr.sh input.vcf(.gz)
