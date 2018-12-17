#!/bin/bash

## This script takes a gff3 and a chromSizes file and extracts the sub-features in BED format:
## (exon, intron, intergenic, CDS, UTR, transcript, gene).

if [[ $# -eq 0 ]] ; then
    echo 'Usage: ./GenomicAnnotate_prepare.sh in.gff3 chromSizes.txt'
    exit 0
fi

#################################################################################################################################
#################################################################################################################################

## Sorting the gff3 and chromSizes files:
cat $1 | mawk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n --parallel=16"}' > ${1%.gff3}_sorted.gff3
sort -k1,1 -k2,2n $2 > ${2%.txt}_sorted.txt

#################################################################################################################################
#################################################################################################################################
## Extract the different features of the gff3 into separate files

## Do with this function:
SPLIT_GFF3 () {
  FILE=$1
  mawk -v FEATURE=$2 '$1 ~ /^#/ {print $0;next} {if ($3 == FEATURE) print $0}' $FILE 
}

## Run functions for exon, CDS, UTRs, transcript, gene:
SORTED_GFF=${1%.gff3}_sorted.gff3
SPLIT_GFF3 $SORTED_GFF exon | gff2bed < /dev/stdin > ${1%.gff3}_exon_sorted.bed
SPLIT_GFF3 $SORTED_GFF five_prime_UTR | gff2bed < /dev/stdin > ${1%.gff3}_5UTR_sorted.bed
SPLIT_GFF3 $SORTED_GFF three_prime_UTR | gff2bed < /dev/stdin > ${1%.gff3}_3UTR_sorted.bed
SPLIT_GFF3 $SORTED_GFF CDS | gff2bed < /dev/stdin > ${1%.gff3}_CDS_sorted.bed
SPLIT_GFF3 $SORTED_GFF gene | gff2bed < /dev/stdin > ${1%.gff3}_gene_sorted.bed
SPLIT_GFF3 $SORTED_GFF transcript | gff2bed < /dev/stdin > ${1%.gff3}_transcript_sorted.bed

## Add a file for intergenic regions (intergenic is the complement of the entire genome with "gene" feature),
## write as BED in the same format as the output of gff2bed:
bedtools complement -i ${1%.gff3}_gene_sorted.bed -g ${2%.txt}_sorted.txt | \
  mawk 'OFS="\t" {print $1, $2, $3, ".", ".", ".", ".", "intergenic", ".", "."}' > ${1%.gff3}_intergenic_sorted.bed

## Add file with introns (intron is the complement of cat(intergenic, exon) and the rest of the genome).
## Pipe the output of complement into intersect (intersecting with "gene" file) to get the gene that the intron belongs to:
bedtools complement -i <(cat ${1%.gff3}_exon_sorted.bed ${1%.gff3}_intergenic_sorted.bed | sort -k1,1 -k2,2n) -g ${2%.txt}_sorted.txt | \
  bedtools intersect -wa -wb -a - -b ${1%.gff3}_gene_sorted.bed | \
  mawk 'OFS="\t" {print $1, $2, $3, $7, $8, $9, $10, "intron", "." ,$13 }' > ${1%.gff3}_intron_sorted.bed
  
## Make a BED with the first bp of every gene for later use with bedtools closest,
## of course strand specific:
mawk 'OFS="\t" {if ($6 == "+") $3=$2+1; print $0}' ${1%.gff3}_gene_sorted.bed | mawk 'OFS="\t" {if ($6 == "-") $2=$3-1; print $0}' | \
  sort -k1,1 -k2,2n > ${1%.gff3}_gene_firstBp_sorted.bed
