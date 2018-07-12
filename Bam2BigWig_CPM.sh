#!/bin/bash

######################################################################################################
##
## Take an indexed BAM file and create a CPM-normalized bigwig without intermediate files, using
## - sambamba depth +  mawk to get the bedgraph,
## - mawk for custom to-CPM conversion and
## - bg2bw (git: cancerit/cgpBigWig) to read the bg from stdin and output the bw
##
######################################################################################################

BAM=$1
CHROMSIZES=$2

######################################################################################################

## Check bam file:

if [[ $# -ne 2 ]] ; then
    echo 'Usage: ./Bam2BigWig_CPM.sh in.bam chromSizes.txt'
    exit 1
fi

DO_EXIT=$(echo 'BAM file looks corrputed -- exiting' && exit 1)
samtools quickcheck $1 && echo '[INFO]: Processing' $1 || ${DO_EXIT}

## Check if indexed:
if [[ ! -f ${1}.bai ]]; then
  echo '[INFO]: BAM file is not indexed -- indexing now:'
  sambamba index -t 8 $1
  fi

######################################################################################################

## Functions:
SCALE_FACTOR=$(bc <<< "scale=8;1000000/$(samtools idxstats $1 | awk '{SUM+=$3} END {print SUM}')") 

## Depth with sambamba:
sambamba depth base -t 8 $1 | mawk 'NR > 1, OFS="\t" {print $1, $2, $3}' | \
  
  ## Convert depth to bedgraph (nathanhaigh/depth2bedgraph.awk)
  mawk 'BEGIN{FS="\t";OFS="\t"}{
  if(NR>1 && ($1!=prev_seq || $3!=prev_cov || $2>prev_pos+1)){
    print prev_seq,start,end,prev_cov
    start = $2-1
  } else if(NR==1){
    start = $2-1
  }
  end = $2
  prev_seq = $1
  prev_pos = $2
  prev_cov = $3
  }
  END{
    print prev_seq,start,end,prev_cov
  }' | \
    
    ## Normalize:
    mawk -v SF=${SCALE_FACTOR} 'OFS="\t" {print $1, $2, $3, $4*SF}' | \

      ## write to bw (cancerit/cgpBigWig)
      bg2bw -i /dev/stdin -c $CHROMSIZES -o ${1%.bam}_CPM.bigwig  
    
######################################################################################################              
