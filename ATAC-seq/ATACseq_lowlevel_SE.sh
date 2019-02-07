#!/bin/bash

######################################################################################################################################

## Same as ATAC-seq_lowlevel.sh version a5d9825 but for single-end data

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
  totalReads=$(samtools idxstats $1 | awk '{SUM += $3} END {print SUM}')

  echo '[mtDNA Content]' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%' > ${1%.bam}_mtDNA.txt
}; export -f mtDNA

######################################################################################################################################

function Fq2Bam {
  
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
  xargs sambamba view -h -f sam -t 2 --num-filter=0/2308 --filter='mapping_quality > 19' \
    -o /dev/stdout ${BASENAME}_raw.bam | \
    tee >(sam2bed --do-not-sort -R < /dev/stdin | \
          mawk 'OFS="\t" {if ($6 == "+") print $1, $2, $2+1, $4, $5, $6} {if ($6 == "-") print $1, $3-1, $3, $4, $5, $6}' | \
          bedtools slop -s -l 0 -r 159 -g tmp_chromSizes.txt | bgzip -@ 4 > ${BASENAME}_dup.bed.gz) | \
  sambamba view -S -f bam -l 5 -t 6 -o /dev/stdout /dev/stdin | \
    tee ${BASENAME}_dup.bam | \
  sambamba view -l 5 -f bam -t 8 --num-filter=/256 -o ${BASENAME}_dedup.bam /dev/stdin
  
  ls *dup.bam | parallel "sambamba flagstat -t 8 {} > {.}.flagstat"
  
  BamCheck ${BASENAME}_dup.bam
  
  BamCheck ${BASENAME}_dedup.bam
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'ended on') <(date))
  
}; export -f Fq2Bam  
  
####################################################################################################################################

## fastqc:
ls *fastq.gz | parallel "fastqc -t 2 {}"

## Run pipeline:
ls *.fastq.gz | awk -F ".fastq.gz" '{print $1}' | parallel -j 4 "Fq2Bam {} mm10 2>> {}.log"

## Get tracks, unscaled, will later be adjusted with DESeq2 size factors:
ls *_dedup.bam | awk -F "_dedup.bam" '{print $1}' | \
 parallel -j 4 "bamCoverage --bam {}_dedup.bam -o {}_dedup.bam_unscaled.bigwig -bs 1 -p 16 -e 160 2>> {}.log"

## Library Complexity for the full dataset without subsetting to peak regions:
ls *_dup.bed.gz | awk -F ".bed" '{print $1}' | \
 parallel "bgzip -c -d -@ 8 {}.bed.gz | sort -S8G -k1,1 -k2,2n --parallel=8 | preseq c_curve -s 5e+05 -o {}.ccurve /dev/stdin"

## Summary:
multiqc -o multiqc_all ./
