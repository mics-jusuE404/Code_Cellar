#!/bin/bash

## Download SRAs based on a provided list of "\n"-separated IDs (SRRXXXXX) and convert to fastq
## Script assumes to be in the directory where the *.sra are loaded to (as defined in vdb-config)
## Requires prefetch and parallel-fastq-dump (from R.Valiers Github) in $PATH

#########################################################################################################

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=71
#SBATCH --partition=normal
#SBATCH --time=160:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=sra2fq
#SBATCH --output=sra2fq.log

#########################################################################################################

function LoadDump {
  prefetch -X 999999999 $1 
  
  if [[ -e ${1}.sra ]]; then
    parallel-fastq-dump -s ${1}.sra -t 10 -O ./ --tmpdir ./ --split-3 --gzip && rm ${1}.sra
  else
    echo '[ERROR]' $1 'apparently not successfully loaded' && exit 1
  fi
}; export -f LoadDump

cat IDs.txt | parallel -j 7 "LoadDump {}"
