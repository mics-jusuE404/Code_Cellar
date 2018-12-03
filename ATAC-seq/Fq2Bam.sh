#!/bin/bash

######################################################################################################################################

## ATAC-seq - alignment of fastq files to hg38, including adapter trimming without intermediate files
## Assumes script in same dir as fastqs, bwa - samtools - samblaster - cutadapt in PATH

######################################################################################################################################

#######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=70
#SBATCH --partition=hims
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=ATACseq_Align
#######

######################################################################################################################################
######################################################################################################################################

function ExitBam {

  (>&2 echo '[ERROR]' $1 'looks suspicious or is empty -- exiting') && exit 1
  
}; export -f ExitBam

######################################################################################################################################
######################################################################################################################################

## Check if BAM is corrupted or incomplete:
function BamCheck {
  
  BASENAME=${1%_*}
  samtools quickcheck -q $1 && echo '' >/dev/null || ExitBam $1
  
  ## Also check if file is not empty:
  if [[ $(samtools view $1 | head -n 1 | wc -l) < 1 ]]; then
    ExitBam $BASENAME
    fi
  
}; export -f BamCheck  

######################################################################################################################################
######################################################################################################################################

function mtDNA {

  if [[ ! -e ${BASENAME}_raw.bam ]];
      then
      echo '[ERROR] Input BAM for mtDNA content is missing -- exit 1' && exit 1
      fi
  
  if [[ ! -e ${BASENAME}_raw.bam.bai ]];
    then
    sambamba index ${BASENAME}_raw.bam
    fi
    
  mtReads=$(samtools idxstats ${BASENAME}_raw.bam | grep 'chrM' | cut -f 3)
  totalReads=$(samtools idxstats ${BASENAME}_raw.bam | awk '{SUM += $3} END {print SUM}')

  echo '[mtDNA Content]' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%' > ${BASENAME}_mtDNA.txt
}; export -f mtDNA

######################################################################################################################################
######################################################################################################################################

function Fq2Bam {

  BASENAME=$1
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'started on') <(date))
  
  if [[ ! -e ${BASENAME}_1.fastq.gz ]] || [[ ! -e ${BASENAME}_2.fastq.gz ]]; then
    echo '[ERROR] At least one input files is missing -- exiting' && exit 1
    fi
  
  ## Nextera adapter:
  ADAPTER1="CTGTCTCTTATACACATCT"
  ADAPTER2="CTGTCTCTTATACACATCT"
  
  BWA_IDX=/scratch/tmp/a_toen03/Genomes/hg38/bwa_index_noALT_withDecoy/hg38_noALT_withDecoy.fa
  
  
  
  seqtk mergepe ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz | \
    cutadapt -j 4 -a $ADAPTER1 -A $ADAPTER2 --interleaved -m 18 --max-n 0.1 - | \
    bwa mem -v 2 -R '@RG\tID:'${BASENAME}'_ID\tSM:'${BASENAME}'_SM\tPL:Illumina' -p -t 16 ${BWA_IDX} /dev/stdin | \
    samtools fixmate -m -@ 2 -O SAM - - | \
    samblaster --ignoreUnmated 2>> ${BASENAME}.log | \
    sambamba view -f bam -S -l 0 -t 4 -o /dev/stdout /dev/stdin | \
      tee >(sambamba flagstat -t 2 /dev/stdin > ${BASENAME}_raw.flagstat) | \
    sambamba sort -m 4G --tmpdir=./ -l 6 -t 16 -o ${BASENAME}_raw.bam /dev/stdin  
    
    BamCheck ${BASENAME}_raw.bam
    mtDNA ${BASENAME}_raw.bam
    
    samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -vE 'chrM|_random|chrU|chrEBV|\*' | \
      xargs sambamba view -f bam -t 8 --num-filter=1/1284 --filter='mapping_quality > 19' \
        -o /dev/stdout ${BASENAME}_raw.bam | \
        tee ${BASENAME}_sorted.bam | \
      sambamba flagstat -t 4 /dev/stdin > ${BASENAME}_sorted.flagstat 
     
  BamCheck ${BASENAME}_sorted.bam
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'ended on') <(date))
}

export -f Fq2Bam

####################################################################################################################################
####################################################################################################################################

## Run:
ls *_1.fastq.gz | awk -F "_1.fastq.gz" '{print $1}' | parallel -j 4 "Fq2Bam {} 2>> {}.log"
  
## Bigwigs:
ls *_sorted.bam | parallel -k -j 8 "bamCoverage -e --normalizeUsing CPM -bs 1 --bam {} -o {.}_CPM.bigwig -p 16"  
