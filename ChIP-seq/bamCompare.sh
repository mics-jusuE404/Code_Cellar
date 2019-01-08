#!/bin/bash

## Take a replicate of ChIP-seq BAM files and use bamCOmpare from deeptools to create a
## bigwig with the mean signal normalized by SES:

##################
function bamComp {

  bamCompare \
    --bamfile1 ${1}*BiolRep1*_sorted.bam \
    --bamfile2 ${1}*BiolRep2*_sorted.bam \
    --outFileName ${1}_SES.bigwig \
    --outFileFormat bigwig \
    --scaleFactorsMethod SES \
    --operation mean \
    --pseudocount 1 \
    --binSize 10 \
    --blackListFileName /scratch/tmp/a_toen03/Genomes/hg38/Blacklists/hg38.blacklist.bed \
    --numberOfProcessors 36 \
    --extendReads 400
 }; export -f bamComp

bamComp $BASENAME

ls ../$i/*_sorted.bam | awk -F "_BiolRep" '{print $1}' | sort -u | parallel -j 2 "bamComp {}"
