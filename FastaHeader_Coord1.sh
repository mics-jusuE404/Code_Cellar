#!/bin/bash

## Takes a fasta file, created with BEDtools getfasta from a BED file (so header has 0-based coords,
## and outputs the same file with start coordinate + 1 to be compatible with FIMO --parse-genomic-coord,
## which assumes 1-based coordinates in the header

if [[ $# -eq 0 ]] ; then
    echo '[ERROR]: No input file!'
    exit 0
fi

awk -F ":" '{OFS=""; split($2,a,"-"); if(a[1]) print $1,":",a[1]+1,"-",a[2]; else print; }' $1
