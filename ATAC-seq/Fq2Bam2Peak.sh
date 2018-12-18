#!/bin/bash

######################################################################################################################################

## Script for lowlevel processing of ATAC-seq data, assuming paired-end data with names *_1.fastq.gz and *_2.fastq.gz

######################################################################################################################################

#######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=70
#SBATCH --partition=hims
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=ATACseq_Fq2Bam
#######

######################################################################################################################################
######################################################################################################################################

## Exit function if BAM file looks corrupted or is missing:
function ExitBam {

  (>&2 echo '[ERROR]' $1 'looks suspicious or is empty -- exiting') && exit 1
  
}; export -f ExitBam

######################################################################################################################################
######################################################################################################################################

## Check if BAM is corrupted or empty:
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

## Function to get the % of reads mapped to chrM:
function mtDNA {

  
  BamCheck $1
  
  if [[ ! -e ${1}.bai ]];
    then
    sambamba index -t 8 ${1}
    fi
    
  mtReads=$(samtools idxstats $1 | grep 'chrM' | cut -f 3)
  totalReads=$(samtools idxstats $1 | awk '{SUM += $3} END {print SUM}')

  echo '[mtDNA Content]' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%' > ${1%.bam}_mtDNA.txt
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
    samblaster --ignoreUnmated | \
    sambamba view -f bam -S -l 0 -t 4 -o /dev/stdout /dev/stdin | \
      tee >(sambamba flagstat -t 2 /dev/stdin > ${BASENAME}_raw.flagstat) | \
    sambamba sort -m 4G --tmpdir=./ -l 6 -t 16 -o ${BASENAME}_raw.bam /dev/stdin  
    
    BamCheck ${BASENAME}_raw.bam
    mtDNA ${BASENAME}_raw.bam
    
    ## Remove non-primary chromosomes and duplicates:
    samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -vE 'chrM|_random|chrU|chrEBV|\*' | \
      xargs sambamba view -f bam -t 8 --num-filter=1/1284 --filter='mapping_quality > 19' \
        -o /dev/stdout ${BASENAME}_raw.bam | \
        tee ${BASENAME}_sorted.bam | \
        tee >(sambamba flagstat -t 2 /dev/stdin > ${BASENAME}_sorted.flagstat) | \
      sambamba index -t 4 /dev/stdin ${BASENAME}_sorted.bam.bai
     
  BamCheck ${BASENAME}_sorted.bam
  
  ## Browser track only for visualization, therefore RPM-normalized instead of using the DESeq2 scaling factor:
  bamCoverage -bs 1 -p 16 -e --normalizeUsing CPM -o ${BASENAME}_sorted_CPM.bigwig --bam ${BASENAME}_sorted.bam
  
  ## Clean up:
  if [[ ! -d BAM_raw ]]; then
    mkdir BAM_raw; fi
    mv ${BASENAME}_raw* BAM_raw
  
  if [[ ! -d BAM_sorted ]]; then
    mkdir BAM_sorted; fi
    mv ${BASENAME}_sorted* BAM_raw
  
  if [[ ! -d fastq ]]; then
    mkdir fastq; fi
    mv ${BASENAME}*.fastq.gz fastq
    
  if [[ ! -d logs ]]; then
    mkdir logs; fi
    mv ${BASENAME}*.log logs
    mv {BASENAME}_mtDNA.txt logs
    
  if [[ ! -d bigwig ]]; then
    mkdir bigwig; fi
    mv ${BASENAME}*bigwig bigwig
    
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'ended on') <(date))
  
}

export -f Fq2Bam

####################################################################################################################################
####################################################################################################################################

## Run:
ls *_1.fastq.gz | awk -F "_1.fastq.gz" '{print $1}' | parallel -j 4 "Fq2Bam {} 2>> {}.log"

## Call peaks with default FDR settings, can be filtered more stringently lateron:
source activate py27
ls *_sorted.bam | \
  awk -F "_sorted.bam" '{print $1}' | \
  parallel "macs2 callpeak -f BAMPE --nomodel --keep-dup=all -g hs -n {} -t {}_sorted.bam"
