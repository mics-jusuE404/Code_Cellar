#!/bin/bash

## Lowlevel workflow for single-end MPRA assuming Nextera-based adapters
## Naming conventions are ${BASENAME}_(R/D)NA_rep*.fastq.gz

#######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=70
#SBATCH --partition=hims
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=STARR_Align
#SBATCH --output=Fq2Bam.log
#######

####################################################################################################################################
####################################################################################################################################

BASENAME=$1

####################################################################################################################################
####################################################################################################################################

## Adapter-trim and align data to hg38:
function Fq2Bam {

  if [[ ! -e ${BASENAME}.fastq.gz ]]; then
    echo '[ERROR] Input file missing -- exiting' && exit 1; fi
  
  ## Adapters/Index:
  ADAPTER1="CTGTCTCTTATACACATCT"
  BWA_IDX=/scratch/tmp/a_toen03/Genomes/hg38/bwa_index_noALT_withDecoy/hg38_noALT_withDecoy.fa
  
  ####################################################################################################################################
  
  echo '[START]' $BASENAME 'on:' >> ${BASENAME}.log && date >> ${BASENAME}.log && echo '' >> ${BASENAME}.log
  
  cutadapt -j 4 -a $ADAPTER1 -m 36 --max-n 0.1 ${BASENAME}.fastq.gz 2>> ${BASENAME}.log | \
    bwa mem -v 2 -R '@RG\tID:'${BASENAME}'_ID\tSM:'${BASENAME}'_SM\tPL:Illumina' -p -t 16 ${BWA_IDX} /dev/stdin 2>> ${BASENAME}.log | \
    samblaster --ignoreUnmated 2>> ${BASENAME}.log | \
    sambamba view -f bam -S -l 1 -t 4 -o /dev/stdout /dev/stdin 2>> ${BASENAME}.log | \
    sambamba sort -m 2G --tmpdir=./ -l 6 -t 16 -o ${BASENAME}_raw.bam /dev/stdin 2>> ${BASENAME}.log
        
    samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -vE 'chrM|_random|chrU|chrEBV|\*' 2>> ${BASENAME}.log | \
      xargs sambamba view -f bam -t 8 --num-filter=0/4 --filter='mapping_quality > 19' 2>> ${BASENAME}.log \
      -o ${BASENAME}_sortedDup.bam ${BASENAME}_raw.bam 2>> ${BASENAME}.log 
    
    sambamba view -f bam -t 8 --num-filter=/1028 -o ${BASENAME}_sortedDeDup.bam ${BASENAME}_sortedDup.bam 2>> ${BASENAME}.log
    
    ls ${BASENAME}*.bam | parallel "sambamba flagstat -t 8 {} > {.}.flagstat 2>> ${BASENAME}.log"
    
  echo '[END]' $BASENAME 'on:' >> ${BASENAME}.log && date >> ${BASENAME}.log
  
}; export -f Fq2Bam

####################################################################################################################################
####################################################################################################################################

## Get the percentage of mitochondrial DNA in library:
function mtDNA {

 mtReads=$(samtools idxstats ${BASENAME}_raw.bam | grep 'chrM' | cut -f 3)
 totalReads=$(samtools idxstats ${BASENAME}_raw.bam | awk '{SUM += $3} END {print SUM}')

 echo '[mtDNA Content]' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%' > ${BASENAME}_mtDNA.txt
 
}; export -f mtDNA

####################################################################################################################################
####################################################################################################################################

## Combine all DNA/RNA files including duplicates into one file to run with preseq:
function COMPLEXITY {
  
  ## Combine using sambamba merge:
  find ./ -maxdepth 1 -name "*_rep*_sortedDup.bam" | \
    awk -F "_rep" '{print $1}' | \
    sort -k1,1 -u | \
    parallel -j 4 "sambamba merge -t 16 {}_combined_sortedDup.bam {}*rep*_sortedDup.bam"
    
  ## Run preseq c_curve:
  find ./ -maxdepth 1 -name "*_combined_sortedDup.bam" | \
    parallel "preseq c_curve -o {.}_ccurve -s 5e+05 -seed 1 {}"
    
} export -f COMPLEXITY

####################################################################################################################################
####################################################################################################################################

## Alignment:
ls *.fastq.gz | awk -F ".fastq.gz" '{print $1}' | \
  parallel -j 4 "Fq2Bam {} && mtDNA {}"
  
## Preseq
COMPLEXITY  
  
## Bigwigs:
ls *_sortedDup.bam | parallel -j 4 "bamCoverage -e --normalizeUsing CPM -bs 1 --bam {} -o {.}_CPM.bigwig -p 16" 
