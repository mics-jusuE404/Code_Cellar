#!/bin/bash

## Extract cutting sites from ATAC-seq data (=the start and end coordinate of every fragment,
## shited +4 if on top strand and -5 if on bottom strand,
## normalize by total number of fragments and write as bigwig:

SCALE=$(bc <<< "scale=6;1000000/$(sambamba view -c $1)")
sambamba sort -n --tmpdir=./ -t 32 -m 40G -o ${1%.bam}_nsort.bam $1 

samtools view -b -@ 8 -L primary.bed ${1%.bam_nsort.bam | \
  bedtools bamtobed -bedpe -i - | \
  mawk 'OFS="\t" { if ($9=="+") {print $1, $2+4, $2+5 "\n" $1, $6+3, $6+4} else if ($9=="-") {print $1, $2-5, $2-4 "\n" $1, $6-6, $6-5}}' | \
  sort -k1,1 -k2,2n --parallel=32 | \
  bedtools genomecov -bg -scale $SCALE -i - -g /scratch/tmp/a_toen03/Genomes/mm10/mm10_chromSizes.txt > ${1%.bam}_cutSites.bedGraph

bedGraphToBigWig ${1%.bam}_cutSites.bedGraph /scratch/tmp/a_toen03/Genomes/mm10/mm10_chromSizes.txt ${1%.bam}_cutSites_cpm.bigwig
