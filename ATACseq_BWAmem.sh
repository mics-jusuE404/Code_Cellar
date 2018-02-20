#!/bin/bash

#### Generic ATAC-seq pipeline for low-level processing from fastq to filtered BAM files using BWA mem for alignment
#### Written by Alexander Toenges (a.toenges@uni-muenster.de) 2017
#### Assumes to be in same dir as the fastq files and bwa, samtools, sambamba , GNU parallel and bamCoverage (deeptools) in PATH:
#### ...and yes the filtering after the alignment is probably fartoo bulky, but I prefer to feed BAMs into MACS for peak calling
#### that really only contain the reads to be used by MACS2 (macs2 callpeak -t BAM --keep-dup=all --nomodel -f BAMPE)
#################################################################################################################################
#################################################################################################################################

## Choose the genome and set if fastq file are compressed or not
## ASSEMBLY can be "hg38", "hg19" or "mm10":
ASSEMBLY="hg38"
COMPRESSED="yes"

BWA_HG38="/path/to/bwa_idx"
BWA_HG19="/path/to/bwa_idx"
BWA_MM10="/path/to/bwa_idx"

#################################################################################################################################
#################################################################################################################################

## Pipeline Start
echo '[INFO]: Started'
date
echo ''

## Set species and exit if no valid assembly was specified:
if [[ $ASSEMBLY == "hg38" ]]; then
  TMP_GENOME=${BWA_HG38}
  ECHOOUT="hg38"
  else
    if [[ $ASSEMBLY == "hg19" ]]; then
      TMP_GENOME=${BWA_HG19}
      ECHOOUT="hg19"
      else
        if [[ $ASSEMBLY == "MM10" ]]; then
          TMP_GENOME=${BWA_MM10}
          ECHOOUT="mm10"
        else
          echo '[ERROR]: Neither hg19/hg38/mm10 selected!'
          exit
        fi  
      fi
    fi  

echo '=> Alignment to' $ECHOOUT
  
# Start pipeline:
for i in *_1.fastq* 
  do 
  BASENAME=${i%_1.fastq*}
  echo ''
  echo '### Currently in process:' $BASENAME
  echo ''
    
  ## Adapter trimming with skewer:
  
  # --- -y and -x are the Nextera transposase adapter sequences
  # --- remove degenerated (NNNNN) reads ::: -n
  # --- trim 3' until hitting a base with qual >= 30 ::: -q 30
  # --- remove reads with mean quality below 25 ::: -Q 25
  # --- multithreaded ::: -t INT
  
  echo '=> Adapter trimming by skewer'
    if [[ $COMPRESSED == "yes" ]]; then
      skewer -x CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -y CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
      -m pe -t 4 --quiet ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz --output ${BASENAME} -n -q 30 -Q 25
      echo ''
      else
        skewer -x CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -y CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
        -m pe -t 8 --quiet ${BASENAME}_1.fastq ${BASENAME}_2.fastq --output ${BASENAME} -n -q 30 -Q 25
        echo ''
      fi 
      
  ###################################################################################################
  ###################################################################################################
    
  # ALIGNMENT WITH BWA MEM --- OUTPUTFILE IS are called *_raw.bam
  echo '=> Alignment:' $BASENAME
  bwa mem -v 2 -t 4 ${TMP_GENOME} ${BASENAME}-trimmed-pair1.fastq ${BASENAME}-trimmed-pair2.fastq | \
    sambamba view -S -f bam -l 0 /dev/stdin | \
    sambamba sort -m 6G --tmpdir=./ -l 9 -t 3 -o ${BASENAME}_raw.bam /dev/stdin
    samtools flagstat ${BASENAME}_raw.bam > ${BASENAME}_raw.txt
  rm ${BASENAME}*.fastq
  
  ###################################################################################################
  ###################################################################################################
  
  # Filtering and sorting:
  # Remove: all chromosomes that are not chr1-22,X (so M, Y, all random and U)
  echo '=> Filtering:' $BASENAME
  KEEPCHR='chr1\|chr2\|chr3\|chr4\|chr5\|chr6\|chr7\|chr8\|chr9\|chr10\|chr11\|chr12\|chr13\|chr14\|chr15\|chr16\|chr17\|chr18\|chr19\|chr20\|chr21\|chr22\|chrX'
                         
  samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -w $KEEPCHR | \
    xargs samtools view -bhu -f 2 -F 2304 -q 30 -@ 3 -o ${BASENAME}_fi.bam ${BASENAME}_raw.bam
    samtools index -@ 3 ${BASENAME}_fi.bam
    samtools flagstat -@ 3 ${BASENAME}_fi.bam > ${BASENAME}_fi.txt
    
  # remove duplicates:
  echo '=> Deduplication:' $BASENAME
  sambamba markdup -r -t 4 ${BASENAME}_fi.bam ${BASENAME}_rmdup.bam
  samtools flagstat -@ 3 ${BASENAME}_rmdup.bam > ${BASENAME}_rmdup.txt
  echo ''
  done
  
  ## Browser tracks:
  echo '=> Bigwigs:'
  ls *_rmdup.bam | parallel "bamCoverage --bam {} -bs 1 --normalizeUsingRPKM --numberOfProcessors --outFileName {.}.bigwig"
