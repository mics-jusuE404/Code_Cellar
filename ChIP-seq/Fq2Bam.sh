#!/bin/bash

######################################################################################################################################

## ChIP-seq - alignment of singlg-end fastq files to hg38, including adapter trimming without intermediate files
## Assumes script in same dir as fastqs, bwa - samtools - samblaster - cutadapt in PATH

######################################################################################################################################

#######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=70
#SBATCH --partition=hims
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=ChIPseq_Align
#SBATCH --output=Fq2Bam.log
#######

function Fq2Bam {

  BASENAME=$1
 
  if [[ ! -e ${BASENAME}.fastq.gz ]]; then
    echo '[ERROR] Input file is missing -- exiting' && exit 1
    fi
  
  ## Nextera adapter:
  ADAPTER1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
  
  BWA_IDX=/scratch/tmp/a_toen03/Genomes/hg38/bwa_index_noALT_withDecoy/hg38_noALT_withDecoy.fa
  
  ####################################################################################################################################
  
  echo '#############################################################' >> ${BASENAME}.log
  echo '[START]' $BASENAME 'on:' >> ${BASENAME}.log && date >> ${BASENAME}.log && echo '' >> ${BASENAME}.log
  
  cutadapt -j 4 -a $ADAPTER1 -m 18 --max-n 0.1 ${BASENAME}.fastq.gz 2>> ${BASENAME}.log | \
    bwa mem -v 2 -R '@RG\tID:'${BASENAME}'_ID\tSM:'${BASENAME}'_SM\tPL:Illumina' -t 16 ${BWA_IDX} /dev/stdin 2>> ${BASENAME}.log | \
    samblaster --ignoreUnmated 2>> ${BASENAME}.log 2>> ${BASENAME}.log | \
    sambamba view -f bam -S -l 1 -t 4 -o /dev/stdout /dev/stdin 2>> ${BASENAME}.log | \
    sambamba sort -m 2G --tmpdir=./ -l 6 -t 16 -o ${BASENAME}_raw.bam /dev/stdin 2>> ${BASENAME}.log
        
    samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -v 'chrM|_random|chrU|chrEBV|\*' 2>> ${BASENAME}.log | \
      xargs sambamba view -f bam -t 8 --num-filter=0/1284 --filter='mapping_quality > 19' 2>> ${BASENAME}.log \
      -o ${BASENAME}_sorted.bam ${BASENAME}_raw.bam 2>> ${BASENAME}.log 
    
    ls ${BASENAME}*.bam | parallel "sambamba flagstat -t 8 {} > {.}.flagstat 2>> ${BASENAME}.log"
    
  echo '[END]' $BASENAME 'on:' >> ${BASENAME}.log && date >> ${BASENAME}.log
  echo '#############################################################' >> ${BASENAME}.log
}

export -f Fq2Bam

## Alignment:
ls *_1.fastq.gz | awk -F "_1.fastq.gz" '{print $1}' | \
  parallel -j 4 "Fq2Bam {} && mtDNA {}"
  
## Bigwigs:
ls *_sorted.bam | parallel -j 4 "bamCoverage -e 400 --normalizeUsing CPM -bs 1 --bam {} -o {.}_CPM.bigwig -p 16"
