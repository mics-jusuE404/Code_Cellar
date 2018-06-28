#!/bin/bash

#### Standard ATAC-seq pipeline for low-level processing from fastq to filtered BAM files using BWA mem for alignment
#### Written by Alexander Toenges (a.toenges@uni-muenster.de) 2017
#### Assumes to be in same dir as the fastq files and bwa, samtools, sambamba , GNU parallel and bamCoverage (deeptools) in PATH:
#### 
#################################################################################################################################
#################################################################################################################################

## Path to BWA index:
TMP_GENOME=
## COmpressed fastq, "yes" or "no"
COMPRESSED=
#################################################################################################################################
#################################################################################################################################

## Pipeline Start
echo '[INFO]: Started'
date
echo ''

## Set species and exit if no valid assembly was specified:

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
  bwa mem -v 2 -t 24 ${TMP_GENOME} ${BASENAME}-trimmed-pair1.fastq ${BASENAME}-trimmed-pair2.fastq | \
    sambamba view -S -f bam -l 0 -t 4 /dev/stdin | \
    sambamba sort -m 50G --tmpdir=./ -l 5 -t 8 -o ${BASENAME}_raw.bam /dev/stdin
    samtools flagstat ${BASENAME}_raw.bam > ${BASENAME}_raw.txt
  rm ${BASENAME}*.fastq
  
  ###################################################################################################
  ###################################################################################################
  
  # Filtering and sorting:
  # Remove: all chromosomes that are not chr1-22,X (so M, Y, all random and U)
  echo '=> Filtering:' $BASENAME
                        
  samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -v 'chrM' | \
    xargs samtools view -bh -f 2 -F 2304 -q 30 -@ 3 -o ${BASENAME}_fi.bam ${BASENAME}_raw.bam
    sambamba index -t 8 ${BASENAME}_fi.bam
    sambamba flagstat -t 8 ${BASENAME}_fi.bam > ${BASENAME}_fi.txt
  
  ###################################################################################################
  ###################################################################################################
  
  # remove duplicates:
  echo '=> Deduplication:' $BASENAME
  sambamba markdup -r -t 4 ${BASENAME}_fi.bam ${BASENAME}_rmdup.bam
  sambamba flagstat -t 8 ${BASENAME}_rmdup.bam > ${BASENAME}_rmdup.txt
  echo ''
  done
  
  ###################################################################################################
  ###################################################################################################
  
  ## Browser tracks:
  echo '=> Bigwigs:'
  ls *_rmdup.bam | parallel "bamCoverage --bam {} -bs 1 --normalizeUsingRPKM --numberOfProcessors 16 --outFileName {.}.bigwig"
