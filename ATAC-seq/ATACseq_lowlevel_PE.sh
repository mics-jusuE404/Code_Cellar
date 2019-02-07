#!/bin/bash

######################################################################################################################################

## Script for lowlevel processing of ATAC-seq data, assuming paired-end data with names *_1.fastq.gz and *_2.fastq.gz
## Don't forget to set appropriate genome version flag!

######################################################################################################################################

#######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --partition=normal
#SBATCH --mem=80G
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=ATACseq_Fq2Bam
#######

BASENAME=$1
 
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
  sambamba view -f sam -t 8 --num-filter=1/2308 --filter='mapping_quality > 19' ${BASENAME}_rawbackup.bam | \
  bedtools bamtobed -bedpe -i - | \
  mawk 'OFS="\t" {if ($1 !~ /chrM|_random|chrU|chrEBV/) print $1, $2, $6, $7, $8, $9}' | \
  sort -S10G -k1,1 -k2,2n --parallel=16 | \
  bgzip -@ 8 > ${BASENAME}_dup_bedpe.bed.gz
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'ended on') <(date))
  
}

export -f Fq2Bam  
  
####################################################################################################################################
####################################################################################################################################

## fastqc:
ls *fastq.gz | parallel "fastqc -t 2 {}"

## Run pipeline:
ls *_1.fastq.gz | awk -F "_1.fastq.gz" '{print $1}' | parallel -j 4 "Fq2Bam {} mm10 2>> {}.log"

## Get browser tracks, unscaled, will be later adjusted with DESeq2 size factors:
ls *_dedup.bam | parallel -j 4 "bamCoverage --bam {} -o {.}_unscaled.bigwig -bs 1 -p 16 -e"

## Library Complexity:
${BASENAME}_dup_bedpe.bed.gz | awk -F ".bed.gz" '{print $1}' | \
 parallel "bgzip -c -d -@ 8 {}.bed.gz | preseq c_curve 5e+05 -o {}_ccurve.txt /dev/stdin"

## Insert Sizes:
conda activate R
ls *_dedup.bam | parallel "picard CollectInsertSizeMetrics I={} O={.}_InsertSizes.txt H={.}_InsertSizes.pdf quiet=true"
conda deactivate
