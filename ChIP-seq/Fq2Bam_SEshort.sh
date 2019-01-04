#!/bin/bash

######################################################################################################################################

## ChIP-seq - alignment of single-end fastq files to hg38 assuming very short reads from ENCODE:
## Assumes script in same dir as fastqs, bwa - samtools - samblaster - deeptools - fastqc
## Last update: 20.11.18

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
  
  ALN_IDX=/scratch/tmp/a_toen03/Genomes/hg38/bowtie_index_noALT_withDecoy/hg38_noALT_withDecoy
  
  ####################################################################################################################################
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'started on') <(date))
  
  ## Map with bowtie trimming to 36bp:
  seqtk trimfq -L 36 ${BASENAME}.fastq.gz | \
    bowtie --sam --quiet -m 1 --best --strata --threads 16 ${ALN_IDX} /dev/stdin | \
    samblaster --ignoreUnmated | \
      tee >(samtools flagstat /dev/stdin > ${BASENAME}_raw.flagstat) | \
    sambamba view -f bam -S -l 1 -t 4 -o /dev/stdout /dev/stdin | \
    sambamba sort -m 4G --tmpdir=./ -l 6 -t 16 -o ${BASENAME}_raw.bam /dev/stdin
  
  ## Remove unmapped and duplicated reads (1028), also index/flagstat:
  samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -vE 'chrM|_random|chrU|chrEBV|\*' | \
    xargs sambamba view -f bam -t 8 --num-filter=0/1028 --filter='mapping_quality > 19' \
      -o /dev/stdout ${BASENAME}_raw.bam | \
        tee >(samtools index -@ 2 - ${BASENAME}_sorted.bam.bai) | \
        tee ${BASENAME}_sorted.bam | \
      samtools flagstat -@ 2 /dev/stdin > ${BASENAME}_sorted.flagstat
        
    
  ls ${BASENAME}*.bam | parallel "sambamba flagstat -t 8 {} > {.}.flagstat"
  
  bamCoverage --normalizeUsing CPM --bam ${BASENAME}_sorted.bam -o ${BASENAME}_sorted_CPM.bigwig -bs 1 -e 500
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'ended on') <(date))
    
}

export -f Fq2Bam

## Alignment:
ls *.fastq.gz | awk -F ".fastq.gz" '{print $1}' | parallel -j 4 "Fq2Bam {} 2>> {}.log"
