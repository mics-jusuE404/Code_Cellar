#!/bin/bash

## Extract cutting sites from ATAC-seq data (=the start and end coordinate of every fragment,
## shited +5 if on top strand and -4 if on bottom strand,
## normalize by total number of fragments and write as bigwig:

SCALE=$(bc <<< "scale=6;1000000/$(sambamba view --num-filter=64/ -c $1)")
samtools view -f 64 -@ 4 $1 | \
  mawk 'OFS="\t" { if ($9 > 0) {print $3, $4-1+5, $4-1+5+1 "\n" $3, $4-1+5+$9, $4-1+5+$9+1} \
    else if ($9 < 0) {print $3, $4-1-4, $4-1-4+1, "\n" $3, $4-1-4+$9, $4-1-4+$9+1}}' | \
      bedtools genomecov -i - -bg -scale $SCALE -g chromSizes.txt > ${1%.bam}_cutsites_cpm.bg
        
bedGraphToBigWig ${1%.bam}_cutsites_cpm.bg chromSizes.txt ${1%.bam}_cutsites_cpm.bigwig
