#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --partition=normal
#SBATCH --mem=80G
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=ATACseq_Fq2Bam

######################################################################################################################################

GENOME="mm10"
MODE="PE"

######################################################################################################################################

## Script for lowlevel processing of ATAC-seq data:
## Don't forget to set appropriate genome version flag!

######################################################################################################################################
## print script for sizefactors:

echo '
## Script to calculate DESeq2 size factors for ATAC-seq data based on a count matrix
## over the entire genome with 500bp windows, taking the top 100k windows based on rowMeans.
## File comes from stdin via the main bash script ATACseq_lowlevel(...).sh

suppressMessages(library(DESeq2))
suppressMessages(library(data.table))
 
## read data:
counts <- fread('cat /dev/stdin', skip = 1, header = T, data.table = F)

## get rowMeans
counts <- data.frame(counts, rowMean = rowMeans(counts[,7:ncol(counts)]))

## sort by rowMeans
counts <- counts[with(counts, order(-rowMean)),]

## top 100k regions based on rowMeans:
counts_subset <- head(counts[,7:(ncol(counts)-1)], n=100000)

## size factors:
SizeFactors <- estimateSizeFactorsForMatrix(counts = counts_subset)

## write:
write.table(data.frame(SizeFactors), sep="\t", col.names = F, row.names = T,
quote = F, file = "./sizeFactors.txt")
' > sizeFactors.R && chmod +x sizeFactors.R

######################################################################################################################################

## Exit function if BAM file looks corrupted or is missing:
function ExitBam {

  (>&2 echo '[ERROR]' $1 'looks suspicious or is empty -- exiting') && exit 1
  
}; export -f ExitBam

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

## Function to get the % of reads mapped to chrM:
function mtDNA {

  
  BamCheck $1
  
  if [[ ! -e ${1}.bai ]];
    then
    sambamba index -t 8 ${1}
    fi
    
  mtReads=$(samtools idxstats $1 | grep 'chrM' | cut -f 3)
  totalReads=$(samtools idxstats $1 | awk '{SUM += $MODE} END {print SUM}')

  echo '[mtDNA Content]' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%' > ${1%.bam}_mtDNA.txt
}; export -f mtDNA

######################################################################################################################################

function Fq2BamPE {

 BASENAME=$1

 (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'started on') <(date))

 if [[ ! -e ${BASENAME}_1.fastq.gz ]] || [[ ! -e ${BASENAME}_2.fastq.gz ]]; then
   echo '[ERROR] At least one input files is missing -- exiting' && exit 1
   fi

 ## Nextera adapter:
 ADAPTER1="CTGTCTCTTATACACATCT"
 ADAPTER2="CTGTCTCTTATACACATCT"

 if [[ ${2} == "hg38" ]]; then
   BWA_IDX="/scratch/tmp/a_toen03/Genomes/hg38/bwa_index_noALT_withDecoy/hg38_noALT_withDecoy.fa"
   fi

 if [[ ${2} == "mm10" ]]; then
   BWA_IDX="/scratch/tmp/a_toen03/Genomes/mm10/bwa_index/mm10.fa"
   fi  

 seqtk mergepe ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz | \
 cutadapt -j 4 -a $ADAPTER1 -A $ADAPTER2 --interleaved -m 18 --max-n 0.1 - | \
 bwa mem -v 2 -R '@RG\tID:'${BASENAME}'_ID\tSM:'${BASENAME}'_SM\tPL:Illumina' -p -t 16 ${BWA_IDX} /dev/stdin | \
 samtools fixmate -m -@ 2 -O SAM - - | \
 samblaster --ignoreUnmated | \
 sambamba view -f bam -S -l 0 -t 4 -o /dev/stdout /dev/stdin | \
   tee ${BASENAME}_rawbackup.bam | \
   tee >(sambamba flagstat -t 2 /dev/stdin > ${BASENAME}_raw.flagstat) | \
 sambamba sort -m 4G --tmpdir=./ -l 6 -t 16 -o ${BASENAME}_raw.bam /dev/stdin  

 BamCheck ${BASENAME}_raw.bam
 mtDNA ${BASENAME}_raw.bam

 ## Remove non-primary chromosomes and duplicates:
 samtools idxstats ${BASENAME}_raw.bam | tee tmp_chromSizes.txt | cut -f 1 | grep -vE 'chrM|_random|chrU|chrEBV|\*' | \
  xargs sambamba view -l 5 -f bam -t 8 --num-filter=1/2308 --filter='mapping_quality > 19' \
   -o /dev/stdout ${BASENAME}_raw.bam | \
  tee ${BASENAME}_dup.bam | \
  tee >(samtools index - ${BASENAME}_dup.bam.bai) | \
  sambamba view -l 5 -f bam -t 8 --num-filter=/256 -o ${BASENAME}_dedup.bam /dev/stdin

 ls *dup.bam | parallel "sambamba flagstat -t 8 {} > {.}.flagstat"

 BamCheck ${BASENAME}_dedup.bam

 ## Bedpe:
 sambamba view -f bam -l 0 -t 8 --num-filter=1/2308 --filter='mapping_quality > 19' ${BASENAME}_rawbackup.bam | \
 bedtools bamtobed -bedpe -i - 2> /dev/null | \
 mawk 'OFS="\t" {if ($1 !~ /chrM|_random|chrU|chrEBV/) print $1, $2, $6, $7, $8, $9}' | \
 sort -S10G -k1,1 -k2,2n --parallel=16 | \
 bgzip -@ 8 > ${BASENAME}_dup.bed.gz
 tabix -p bed ${BASENAME}_dup.bed.gz

 (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'ended on') <(date))

 }

 export -f Fq2BamPE

######################################################################################################################################

function Fq2BamSE {
  
BASENAME=$1
(>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'started on') <(date))

if [[ ! -e ${BASENAME}.fastq.gz ]] ; then
  echo '[ERROR] At least one input files is missing -- exiting' && exit 1
  fi
  
## Nextera adapter:
ADAPTER1="CTGTCTCTTATACACATCT"
    
if [[ ${2} == "hg38" ]]; then
  BWA_IDX="/scratch/tmp/a_toen03/Genomes/hg38/bwa_index_noALT_withDecoy/hg38_noALT_withDecoy.fa"
  fi
  
if [[ ${2} == "mm10" ]]; then
  BWA_IDX="/scratch/tmp/a_toen03/Genomes/mm10/bwa_index/mm10.fa"
  fi  
  
## Alignment, save BAM with all reads included (_raw.bam)
cutadapt -j 4 -a $ADAPTER1 -m 18 --max-n 0.1 ${BASENAME}.fastq.gz | \
bwa mem -v 2 -R '@RG\tID:'${BASENAME}'_ID\tSM:'${BASENAME}'_SM\tPL:Illumina' -t 16 ${BWA_IDX} /dev/stdin | \
samblaster --ignoreUnmated | \
sambamba view -f bam -S -l 0 -t 4 -o /dev/stdout /dev/stdin | \
  tee >(sambamba flagstat -t 2 /dev/stdin > ${BASENAME}_raw.flagstat) | \
sambamba sort -m 4G --tmpdir=./ -l 6 -t 16 -o ${BASENAME}_raw.bam /dev/stdin  
  
BamCheck ${BASENAME}_raw.bam
mtDNA ${BASENAME}_raw.bam
  
## Remove non-primary chromosomes and low-quality alignments, 
## outputting BAMs with and w/o duplicates, 
## also output in BED format elongated to 160bp for average fragment size
samtools idxstats ${BASENAME}_raw.bam | tee tmp_chromSizes.txt | cut -f 1 | grep -vE 'chrM|_random|chrU|chrEBV|\*' | \
xargs sambamba view -f bam -l 5 -t 8 --num-filter=0/2308 --filter='mapping_quality > 19' \
  -o /dev/stdout ${BASENAME}_raw.bam | \
tee ${BASENAME}_dup.bam | \
sambamba view -l 5 -f bam -t 8 --num-filter=/256 -o ${BASENAME}_dedup.bam /dev/stdin
  
ls *dup.bam | parallel "sambamba flagstat -t 8 {} > {.}.flagstat"
  
## BED file for dup reads to be used with preseq:
bedtools bamtobed -i ${BASENAME}_dup.bam | \
mawk 'OFS="\t" {if ($6 == "+") print $1, $2, $2+1, $4, $5, $6} {if ($6 == "-") print $1, $MODE-1, $MODE, $4, $5, $6}' | \
bedtools slop -s -l 0 -r 159 -g tmp_chromSizes.txt | \
sort -S10G -k1,1 -k2,2n --parallel=16 | \
bgzip -@ 4 > ${BASENAME}_dup.bed.gz
tabix -p bed ${BASENAME}_dup.bed.gz
    
BamCheck ${BASENAME}_dup.bam
  
BamCheck ${BASENAME}_dedup.bam
  
(>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'ended on') <(date))
  
}; export -f Fq2BamSE

####################################################################################################################################

RSCRIPT="$HOME/anaconda3_things/anaconda3/envs/R_env/bin/Rscript"
function SizeFactor {

 ## 500bp windows over the genome:
 bedtools makewindows -w 500 -g tmp_chromSizes.txt | \
 mawk 'OFS="\t" {print $1$2$3, $1, $2+1, $3, "+"}' > genome_windows.saf
   
 ## count matrix
 featureCounts --read2pos 5 -a genome_windows.saf -F SAF -T 8 -o genome_windows_counts.txt *_dedup.bam
 
 cat genome_windows_counts.txt | $RSCRIPT sizeFactors.R
 
 }; export -f SizeFactor
   
####################################################################################################################################

function Bigwig {

 FILE=$1
 MODUS=$2
 
 ## bamCoverage multiplies with the size factor but deseq2 divides so invert the deseq factor:
 FACTOR=$(grep $FILE sizeFactors.txt | cut -f2 | bc <<< "scale=6;$(grep $FILE sizeFactors.txt | cut -f2)^-1")
 
 if [[ $MODUS == "PE" ]]; then
  bamCoverage --bam $FILE --scaleFactor $FACTOR -e -o ${FILE%.bam}_geoMean.bigwig -p 16 -bs 1
  fi
  
 if [[ $MODUS == "SE" ]]; then
  bamCoverage --bam $FILE --scaleFactor $FACTOR -e 160 -o ${FILE%.bam}_geoMean.bigwig -p 16 -bs 1
  fi
  
 }; export -f Bigwig  

####################################################################################################################################

## fastqc:
ls *fastq.gz | parallel -j 70 "fastqc -t 1 {}"

####################################################################################################################################

## Run pipeline:
if [[ $MODE != "PE" ]] && [[ $MODE != "SE" ]]; then
 echo '[ERROR] Missing SE/PE information in Fq2Bam function -- exiting'
 exit 1; fi
 
if [[ $MODE == "PE" ]]; then
 ls *_1.fastq.gz | awk -F "_1.fastq.gz" '{print $1}' | parallel -j 4 "Fq2BamPE {} $GENOME 2>> {}.log"
 fi

if [[ $MODE == "SE" ]]; then
 ls *.fastq.gz | awk -F ".fastq.gz" '{print $1}' | parallel -j 4 "Fq2BamSE {} $GENOME 2>> {}.log"
 fi

####################################################################################################################################

## Estimate size factors:
SizeFactor

## Get browser tracks, scaled by the size factor from deseq:
ls *dedup.bam | parallel -j 4 "Bigwig {} $MODE"

####################################################################################################################################

if [[ $MODE == "PE" ]]; then
 ## Insert Sizes:
 ls *_dedup.bam | \
  parallel "picard CollectInsertSizeMetrics I={} O={.}_InsertSizes.txt H={.}_InsertSizes.pdf QUIET=true VERBOSITY=ERROR 2> /dev/null"
 fi

####################################################################################################################################

## Library Complexity:
if [[ $MODE == "PE" ]]; then
 ls *_dup.bed.gz | awk -F ".bed.gz" '{print $1}' | \
  parallel "bgzip -c -d -@ 8 {}.bed.gz | preseq c_curve -s 5e+05 -o {}_ccurve.txt /dev/stdin"
 fi
  
if [[ $MODE == "SE" ]]; then
 ls *_dup.bed.gz | awk -F ".bed" '{print $1}' | \
  parallel "zcat {}.bed.gz | preseq c_curve -s 1e+05 -o {}_ccurve.tsv /dev/stdin"
 fi 

####################################################################################################################################

## Peaks for each sample:
MACS="/anaconda3_things/anaconda3/envs/py27_env/bin/macs2"
ls *dedup.bam | \
 parallel "$MACS callpeak -t {} -n {.} -g mm --extsize 150 --shift -75 --nomodel --keep-dup=all -f BAM"


####################################################################################################################################

## Summary report:
multiqc -o multiqc_all ./ 
