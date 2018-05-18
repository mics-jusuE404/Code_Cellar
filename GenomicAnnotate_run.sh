#!/bin/bash

## Script takes a BED file as input, finds overlapping features from a gff3 file and outputs the distances to the two closest genes.
## Requires that the script < GenomicAnnotate_prepare.sh > has previously been run on the gff3 file to be used.
## Output is a TSV with the genomic position from the BED, the overlapping features and its corresponding genes, and the distance to the two closest gene
## USAGE: ./GenomicAnnotate_run.sh in.bed basename_gff

## EXAMPLE, when using chr11:47361739-47362018 as input with the gencode.v27 (hg38, human) annotations:
## \ chr	start	end	overlapping_feature	overlapping_gene	closest_gene_name	distance
## \ chr11	47361739	47362018	intron	SPI1	MYBPC3	9038
## \ chr11	47361739	47362018	intron	SPI1	SPI1	16558


BEDFILE=$1
GFF_BASENAME=$2

## Take as $1 a BED-like file with regions to find the closest gene:

## 1) Intersect with exon, intron, intergenic, CDS, UTRs:
bedtools intersect -a <(cut -f1-3 $1 | sort -k1,1 -k2,2n -u) -b ${GFF_BASENAME}_exon_sorted.bed \
                                                                ${GFF_BASENAME}_intron_sorted.bed \
                                                                ${GFF_BASENAME}_intergenic_sorted.bed \
                                                                ${GFF_BASENAME}_CDS_sorted.bed \
                                                                ${GFF_BASENAME}_5UTR_sorted.bed \
                                                                ${GFF_BASENAME}_3UTR_sorted.bed \
                                                                -wa -wb -sorted | \
                                                                cut -f 1-3,12,14 > ${BEDFILE}_tmp.intersect.txt    

## Extract the gene name from $5 of the intersect command:
cut -f5 ${BEDFILE}_tmp.intersect.txt | awk -F "gene_name=" '{print $2}' | \
  awk -F ";" '{print $1}' | tr -d '"' > ${BEDFILE}_tmp.intersectGene.txt
  
## add a "." to $4 if feature is intergenic so gene_name field is not blank:  
paste <(cut -f1-4 ${BEDFILE}_tmp.intersect.txt) ${BEDFILE}_tmp.intersectGene.txt | awk 'OFS="\t" {if ($5 == "") $5="."; print}' > ${BEDFILE}_tmp.finalIntersect.txt

## Now use bedtools closest to get the two closest genes:
bedtools closest -k 2 -d -a <(cut -f1-5 ${BEDFILE}_tmp.finalIntersect.txt) -b ${GFF_BASENAME}_gene_firstBp_sorted.bed \
  -sorted > ${BEDFILE}_tmp.closest.txt

## Write final output:
cat <(echo 'chr, start, end, overlapping_feature, overlapping_gene, closest_gene_name, distance' | tr -s ', ' '\t') \
  <(paste \
      <(cut -f1-5 ${BEDFILE}_tmp.closest.txt) \
      <(cut -f15 ${BEDFILE}_tmp.closest.txt | awk -F "gene_name=" '{print $2}' |   awk -F ";" '{print $1}' | tr -d '"') \
      <(cut -f16 ${BEDFILE}_tmp.closest.txt)) | \
        sort -k1,1 -k2,2n -k3,3n -k4,4 -k6,6 -u /dev/stdin | sort -k1,1 -k2,2n -k7,7n
      
## Clean up the tmp files:
rm ${BEDFILE}_tmp.intersect.txt
rm ${BEDFILE}_tmp.intersectGene.txt
rm ${BEDFILE}_tmp.finalIntersect.txt
rm ${BEDFILE}_tmp.closest.txt

exit 0
