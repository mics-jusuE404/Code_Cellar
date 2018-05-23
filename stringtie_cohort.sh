#!/bin/bash

## Use stringtie on an RNA-seq cohort to:
## 1) assembly transcripts (together with reference GTF) from each sample
## 2) merge transcripts into a consensus GTF
## 3) estimate abundance based on the consensus GTF for each sample, creating a count matrix for further analysis.

## Version v_1.0, May 2018

###################################################################################################################
###################################################################################################################

## Assumes that script is in the same directory as the coordinate-sorted RNA-seq BAM files:

if [[ $# -eq 0 ]] ; then
    echo ''
    exit 0
fi
