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

function Fq2Bam {

  BASENAME=$1
 
  if [[ ! -e ${BASENAME}_1.fastq.gz ]] || [[ ! -e ${BASENAME}_2.fastq.gz ]]; then
    echo '[ERROR] At least one input files is missing -- exiting' && exit 1

    fi
  
  ## Nextera adapter:
  ADAPTER1="CTGTCTCTTATACACATCT"
  ADAPTER2="CTGTCTCTTATACACATCT"
  
  BWA_IDX=/scratch/tmp/a_toen03/Genomes/hg38/bwa_index_noALT_withDecoy/hg38_noALT_withDecoy.fa
  
  ####################################################################################################################################
  
  echo '#############################################################' >> ${BASENAME}.log
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'started on') <(date))
  
  seqtk mergepe ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz | \
    cutadapt -j 4 -a $ADAPTER1 -A $ADAPTER2 --interleaved -m 18 --max-n 0.1 - | \
    bwa mem -v 2 -R '@RG\tID:'${BASENAME}'_ID\tSM:'${BASENAME}'_SM\tPL:Illumina' -p -t 16 ${BWA_IDX} /dev/stdin | \
    samtools fixmate -m -@ 2 -O SAM - - | \
    samblaster --ignoreUnmated 2>> ${BASENAME}.log | \
    sambamba view -f bam -S -l 1 -t 4 -o /dev/stdout /dev/stding | \
    sambamba sort -m 4G --tmpdir=./ -l 6 -t 16 -o ${BASENAME}_raw.bam /dev/stdin
        
    samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -vE 'chrM|_random|chrU|chrEBV|\*' | \
      xargs sambamba view -f bam -t 8 --num-filter=1/1284 --filter='mapping_quality > 19' \
      -o ${BASENAME}_sorted.bam ${BASENAME}_raw.bam 
    
    ls ${BASENAME}*.bam | parallel "sambamba flagstat -t 8 {} > {.}.flagstat"
    
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'ended on') <(date))
}

export -f Fq2Bam

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
}

export -f mtDNA

## Alignment:
ls *_1.fastq.gz | awk -F "_1.fastq.gz" '{print $1}' | \
  parallel -j 4 "Fq2Bam {} 2>> {}.log && mtDNA {}"
  
## Bigwigs:
ls *_sorted.bam | parallel -k -j 8 "bamCoverage -e --normalizeUsing CPM -bs 1 --bam {} -o {.}_CPM.bigwig -p 16"  
