#!/bin/bash

## very elaboreate script:

conda activate hicexplorer

## MboI:
findRestSite --fasta mm10.fa --searchPattern GATC -o mm10_MboI.bed
