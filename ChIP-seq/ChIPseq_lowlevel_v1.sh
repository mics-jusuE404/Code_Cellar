#!/bin/bash

######################################################################################################################################

## ChIP-seq - alignment of single-end fastq files:
## Assumes script in same dir as fastqs, bwa - samtools - samblaster - deeptools - fastqc
## Last update: 20.11.18

######################################################################################################################################

GENOME="mm10"
MACS="$HOME/anaconda3_things/anaconda3/envs/py27_env/bin/macs2"

#######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --mem=80G
#SBATCH --partition=normal
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=ChIPseq_Align
#SBATCH --output=Fq2Bam.log
#######

function Fq2Bam {

  BASENAME=$1
  
  if [[ $GENOME == "mm10" ]]; then
    IDX="/scratch/tmp/a_toen03/Genomes/mm10/bowtie2_idx/mm10"
    fi
  
  ####################################################################################################################################
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'started on') <(date))
  
  ## Map with bowtie trimming to 36bp:
  bowtie2 --very-sensitive --threads 16 -x ${IDX} -U ${BASENAME}.fastq.gz | \
  samblaster --ignoreUnmated | \
  sambamba view -f bam -S -l 0 -t 2 -o /dev/stdout /dev/stdin | \
  tee >(sambamba flagstat -t 2 /dev/stdin > ${BASENAME}_raw.flagstat) | \
  sambamba sort -m 4G --tmpdir=./ -l 6 -t 16 -o ${BASENAME}_raw.bam /dev/stdin  
  
  ## Remove unmapped and duplicated reads (1028), also index/flagstat:
  samtools idxstats ${BASENAME}_raw.bam | tee tmp_chromSizes.txt | cut -f 1 | grep -vE 'chrM|_random|chrU|chrEBV|\*' | \
  xargs sambamba view -f bam -l 5 -t 8 --num-filter=0/2308 --filter='mapping_quality > 19' \
  -o /dev/stdout ${BASENAME}_raw.bam | \
  tee >(tee ${BASENAME}_dup.bam | samtools index - ${BASENAME}_dup.bam.bai) | \
  sambamba view -l 5 -f bam -t 8 --num-filter=/1024 -o ${BASENAME}_dedup.bam /dev/stdin
        
  ls ${BASENAME}*dup.bam | parallel "sambamba flagstat -t 8 {} > {.}.flagstat"
  
  bamCoverage --bam ${BASENAME}_dedup.bam -o ${BASENAME}_dedup_unscaled.bigwig -bs 1 -e 400 -p 16
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'ended on') <(date))
    
}

export -f Fq2Bam

## Fastqc:
ls *.fastq.gz | parallel -j 70 "fastqc {}"

## Alignment:
ls *.fastq.gz | awk -F ".fastq.gz" '{print $1}' | parallel -j 4 "Fq2Bam {} hg38 2>> {}.log"

## Macs:
if [[ $GENOME == "mm10" ]]; then GFLAG="mm"; fi
ls *_dedup.bam | parallel "$MACS2 callpeak -t {} -n {.} -g $GFLAG -f BAM"

## Preseq:
(>&2 paste -d " " <(echo '[INFO] LibComplexity started on') <(date))
ls *_dup.bam | parallel "preseq c_curve -bam -s 5e+05 -o {.}_ccurve.txt {}"
(>&2 paste -d " " <(echo '[INFO] LibComplexity ended on') <(date))
