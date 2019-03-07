#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --partition=normal
#SBATCH --mem=80G
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=ATACseq_lowlevel

######################################################################################################################################

## Script for lowlevel processing of ATAC-seq data: Alignment, filtering, quality control, browser tracks:

######################################################################################################################################

## due to different conda envs, some things have to be set explicitely:

GENOME="mm10"
MODE="PE"
BLACKLIST="/scratch/tmp/a_toen03/Genomes/mm10/mm10_consensusBL.bed"
##
RSCRIPT="$HOME/anaconda3_things/anaconda3/envs/R_env/bin/Rscript"
MACS="$HOME/anaconda3_things/anaconda3/envs/py27_env/bin/macs2"

######################################################################################################################################

## Check if required tools are in PATH:
if [[ -e missing_tools.txt ]]; then rm missing_tools.txt; fi

function PathCheck {
  
  if [[ $(command -v $1 | wc -l) == 0 ]]; then 
    echo ${1} >> missing_tools.txt
    fi
  
}; export -f PathCheck

TOOLS=(cat seqtk cutadapt bwa samtools samblaster sambamba tee xargs bedtools mawk bgzip tabix \
       sort paste featureCounts bc bamCoverage parallel fastqc picard preseq multiqc $RSCRIPT \
       bigWigToBedGraph $MACS bg2bw pigz)

for i in $(echo ${TOOLS[*]}); do
  PathCheck $i; done
  
if [[ -e missing_tools.txt ]] && [[ $(cat missing_tools.txt | wc -l | xargs) > 0 ]]; then
  echo '[ERROR] Tools missing in PATH -- see missing_tools.txt' && exit 1
  fi

######################################################################################################################################
## print script for sizefactors:

cat <<EOF > sizeFactors.R
## Script to calculate DESeq2 size factors for ATAC-seq data based on a count matrix
## over the entire genome with 500bp windows, taking the top 100k windows based on rowMeans.
## File comes from stdin via the main bash script ATACseq_lowlevel(...).sh

library(DESeq2)
library(data.table)
 
## read data:
counts <- fread('cat /dev/stdin', skip = 1, header = T, data.table = F)

if (ncol(counts) == 7) {
  write.table(data.frame(c("ONESAMPLE")), sep="\t", col.names = F, row.names = T,
              quote = F, file = "./sizeFactors.txt")
  stop("Only one sample present, so skipping calculation!")              
}

## get rowMeans
counts <- data.frame(counts, rowMean = rowMeans(counts[,7:ncol(counts)]))

## sort by rowMeans
counts <- counts[with(counts, order(-rowMean)),]

## top 100k regions based on rowMeans:
counts_subset <- head(counts[,7:(ncol(counts)-1)], n=100000)

## remove regions with rowMeans beyond the 99.9th percentile, aiming to remove artifact regions without 
## requiring an external blacklist:
counts_subset <- counts_subset[which(rowMeans(counts_subset) <= quantile(rowMeans(counts_subset), .99)),]

## size factors:
SizeFactors <- estimateSizeFactorsForMatrix(counts = counts_subset)

## write:
write.table(data.frame(SizeFactors), sep="\t", col.names = F, row.names = T,
quote = F, file = "./sizeFactors.txt")
EOF

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

 function Fq2BamPE {

  BASENAME=$1
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2BamPE for' $1 'started on') <(date))
  
  if [[ ! -e ${BASENAME}_1.fastq.gz ]] || [[ ! -e ${BASENAME}_2.fastq.gz ]]; then
  echo '[ERROR] At least one input files is missing -- exiting' && exit 1
  fi
  
  ## Nextera adapter:
  ADAPTER1="CTGTCTCTTATACACATCT"
  ADAPTER2="CTGTCTCTTATACACATCT"
  
  if [[ ${2} == "hg38" ]]; then
  IDX="/scratch/tmp/a_toen03/Genomes/hg38/bowtie2_index_noALT_withDecoy/hg38_noALT_withDecoy.fa"
  fi
  
  if [[ ${2} == "mm10" ]]; then
  IDX="/scratch/tmp/a_toen03/Genomes/mm10/bowtie2_idx/mm10"
  fi  
  
  ## Fastq to raw (unfiltered) alignments:
  seqtk mergepe ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz | \
  cutadapt -j 1 -a $ADAPTER1 -A $ADAPTER2 --interleaved -m 18 --max-n 0.1 - | \
  bowtie2 --very-sensitive --threads 16 -x $IDX --interleaved - | \
  samtools fixmate -m -@ 2 -O SAM - - | \
  samblaster --ignoreUnmated | \
  sambamba view -f bam -S -l 5 -t 2 -o /dev/stdout /dev/stdin | \
  tee >(sambamba flagstat -t 2 /dev/stdin > ${BASENAME}_raw.flagstat) | \
  sambamba sort -m 5G --tmpdir=./ -l 5 -t 16 -o ${BASENAME}_raw.bam /dev/stdin  
  
  BamCheck ${BASENAME}_raw.bam
  mtDNA ${BASENAME}_raw.bam
  
  ## 1) remove non-primary chromosomes, low qual. and non-primary alignments, but keep duplicates:
  if [[ ! -e tmp_chromSizes.txt ]]; then
    samtools idxstats ${BASENAME}_raw.bam > tmp_chromSizes.txt
    fi
    
  samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -vE 'chrM|_random|chrU|chrEBV|\*' | \
  xargs sambamba view -l 5 -f bam -t 8 --num-filter=1/2308 --filter='mapping_quality > 19' \
  -o /dev/stdout ${BASENAME}_raw.bam | \
  tee ${BASENAME}_dup.bam | \
  tee >(samtools index - ${BASENAME}_dup.bam.bai) | \
  sambamba view -l 5 -f bam -t 8 --num-filter=/1024 -o /dev/stdout /dev/stdin | \
  tee >(tee ${BASENAME}_dedup.bam | samtools index - ${BASENAME}_dedup.bam.bai) | \
  bedtools bamtobed -i - | \
  mawk 'OFS="\t" {if ($6 == "+") print $1, $2+4, $2+5, ".", ".", $6} {if ($6 == "-") print $1, $3-5, $3-4, ".", ".", $6}' | \
  sort -k1,1 -k2,2n -k3,3n -k6,6 -S10G --parallel=10 |
  tee >(bgzip -@ 6 > ${BASENAME}_cutsites.bed.gz) | \
  bedtools genomecov -bg -i - -g tmp_chromSizes.txt | \
  bg2bw -i /dev/stdin -c tmp_chromSizes.txt -o ${BASENAME}_cutsites_noScale.bigwig
  
  BamCheck ${BASENAME}_dup.bam
  
  ls ${BASENAME}*dup.bam | parallel "sambamba flagstat -t 8 {} > {.}.flagstat"
  
  BamCheck ${BASENAME}_dedup.bam
    
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2BamPE for' $1 'ended on') <(date))

}

export -f Fq2BamPE

######################################################################################################################################

function Fq2BamSE {
  
  BASENAME=$1
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2BamSE for' $1 'started on') <(date))
  
  if [[ ! -e ${BASENAME}.fastq.gz ]] ; then
  echo '[ERROR] At least one input files is missing -- exiting' && exit 1
    fi
    
  ## Nextera adapter:
  ADAPTER1="CTGTCTCTTATACACATCT"
      
  if [[ ${2} == "hg38" ]]; then
  IDX="/scratch/tmp/a_toen03/Genomes/hg38/bowtie2_index_noALT_withDecoy/hg38_noALT_withDecoy.fa"
    fi
  
  if [[ ${2} == "mm10" ]]; then
  IDX="/scratch/tmp/a_toen03/Genomes/mm10/bowtie2_idx/mm10"
    fi  
    
  ## Alignment, save BAM with all reads included (_raw.bam)
  cutadapt -j 1 -a $ADAPTER1 -m 18 --max-n 0.1 ${BASENAME}.fastq.gz | \
  bowtie2 --very-sensitive --threads 16 -x $IDX -U - | \
  samblaster --ignoreUnmated | \
  sambamba view -f bam -S -l 0 -t 2 -o /dev/stdout /dev/stdin | \
  tee >(sambamba flagstat -t 2 /dev/stdin > ${BASENAME}_raw.flagstat) | \
  sambamba sort -m 4G --tmpdir=./ -l 6 -t 16 -o ${BASENAME}_raw.bam /dev/stdin  
    
  BamCheck ${BASENAME}_raw.bam
  mtDNA ${BASENAME}_raw.bam
    
  ## Remove non-primary chromosomes and low-quality alignments, 
  ## outputting BAMs with and w/o duplicates, 

  if [[ ! -e tmp_chromSizes.txt ]]; then
    samtools idxstats ${BASENAME}_raw.bam > tmp_chromSizes.txt
    fi
    
  samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -vE 'chrM|_random|chrU|chrEBV|\*' | \
  xargs sambamba view -l 5 -f bam -t 8 --num-filter=0/2308 --filter='mapping_quality > 19' \
  -o /dev/stdout ${BASENAME}_raw.bam | \
  tee ${BASENAME}_dup.bam | \
  tee >(samtools index - ${BASENAME}_dup.bam.bai) | \
  sambamba view -l 5 -f bam -t 8 --num-filter=/1024 -o /dev/stdout /dev/stdin | \
  tee >(tee ${BASENAME}_dedup.bam | samtools index - ${BASENAME}_dedup.bam.bai) | \
  bedtools bamtobed -i - | \
  mawk 'OFS="\t" {if ($6 == "+") print $1, $2+4, $2+5, ".", ".", $6} {if ($6 == "-") print $1, $3-5, $3-4, ".", ".", $6}' | \
  sort -k1,1 -k2,2n -k3,3n -k6,6 -S10G --parallel=10 |
  tee >(bgzip -@ 6 > ${BASENAME}_cutsites.bed.gz) | \
  bedtools genomecov -bg -i - -g tmp_chromSizes.txt | \
  bg2bw -i /dev/stdin -c tmp_chromSizes.txt -o ${BASENAME}_cutsites_noScale.bigwig
    
  ls ${BASENAME}*dup.bam | parallel "sambamba flagstat -t 8 {} > {.}.flagstat"

  BamCheck ${BASENAME}_dup.bam
    
  BamCheck ${BASENAME}_dedup.bam
    
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2BamSE for' $1 'ended on') <(date))
    
}; export -f Fq2BamSE

######################################################################################################################################

function SizeFactor {
  
  (>&2 paste -d " " <(echo '[INFO] SizeFactor started on') <(date))

  ## 500bp windows over the genome:
  bedtools makewindows -w 500 -g tmp_chromSizes.txt | \
  mawk 'OFS="\t" {print $1"_"$2"_"$3, $1, $2+1, $3, "+"}' > genome_windows.saf
     
  ## count matrix
  featureCounts --read2pos 5 -a genome_windows.saf -F SAF -T 8 -o genome_windows_counts.txt *_dedup.bam
   
  cat genome_windows_counts.txt | $RSCRIPT sizeFactors.R
  
  (>&2 paste -d " " <(echo '[INFO] SizeFactor ended on') <(date))
  
}; export -f SizeFactor
   
####################################################################################################################################

function Bigwig {

  FILE=$1
    
  (>&2 paste -d " " <(echo '[INFO]' 'Browser Tracks for' $1 'started on') <(date))
  
  ## bamCoverage multiplies with the size factor but deseq2 divides so invert the deseq factor:
  if [[ $(cat sizeFactors.txt | wc -l | xargs) > 1 ]]; then
    FACTOR=$(grep ${FILE%_cutsites*}_dedup.bam sizeFactors.txt | \
             cut -f2 | bc <<< "scale=6;$(grep ${FILE%_cutsites*}_dedup.bam sizeFactors.txt | cut -f2)^-1")
    fi
    
  if [[ $(cat sizeFactors.txt | wc -l | xargs) == 1 ]]; then
    FACTOR=1
    fi
  
  ## Bigwig with 80bp fragments to smooth the signal:
  bigWigToBedGraph $FILE /dev/stdout | \
  bedtools slop -b 40 -g tmp_chromSizes.txt -i - | \
  bedtools genomecov -g tmp_chromSizes.txt -i - -bg -scale ${FACTOR} | \
  bg2bw -i /dev/stdin -c tmp_chromSizes.txt -o ${FILE%_noScale.bigwig}_geomean.bigwig
  
  (>&2 paste -d " " <(echo '[INFO]' 'Browser Tracks for' $1 'ended on') <(date))
  
}; export -f Bigwig  

####################################################################################################################################

## Peak calling with Genrich:
function GenR {

  BASENAME=$1
  QVAL=$2
  BLACKLIST=$3
  
  sambamba sort --tmpdir=./ -n -t 8 -m 5G -o /dev/stdout ${BASENAME}_dedup.bam | \
  tee ${BASENAME}_dedup_nsort.bam | \
  Genrich -y -t - -o /dev/stdout -E $BLACKLIST -j -q $QVAL | \
  sort -k1,1 -k2,2n > {}_peaks_Genrich.narrowPeak

}; export -f GenR

####################################################################################################################################

## Function to count reads in a 500bp window around called peak summits:
function FRiP {
  
  PEAKS=$1
  CUTSITES=$2
  
  if [[ -e $1 ]] && [[ -e $2 ]]; then
    TOTAL=$(pigz -c -d -p 4 $CUTSITES | wc -l)
    READ=$(bedtools intersect \
           -a $CUTSITES \
           -b $PEAKS -u \
           -wa -sorted | wc -l)
     paste <(echo ${2%_cutsites.bed.gz}) <(bc <<< "scale=6;$READ/$TOTAL") 
     fi

}; export -f FRiP

####################################################################################################################################
####################################################################################################################################

####################################################################################################################################
####################################################################################################################################

####################################################################################################################################
####################################################################################################################################

## fastqc:
ls *fastq.gz | parallel -j 70 "fastqc -t 1 {}"

####################################################################################################################################

## Run pipeline:
if [[ $MODE != "PE" ]] && [[ $MODE != "SE" ]]; then
  echo '[ERROR] Missing SE/PE information in Fq2Bam function -- exiting'
  exit 1; fi
 
if [[ $MODE == "PE" ]]; then
  ls *_1.fastq.gz | awk -F "_1.fastq.gz" '{print $1}' | sort -u | parallel -j 4 "Fq2BamPE {} $GENOME 2>> {}.log"
  fi

if [[ $MODE == "SE" ]]; then
  ls *.fastq.gz | awk -F ".fastq.gz" '{print $1}' | sort -u | parallel -j 4 "Fq2BamSE {} $GENOME 2>> {}.log"
  fi

####################################################################################################################################

## Estimate size factors:
SizeFactor 2> sizeFactors.log

####################################################################################################################################

## Get browser tracks, scaled by the size factor from deseq:
ls *_cutsites_noScale.bigwig | parallel -j 4 "Bigwig {}"

####################################################################################################################################

## Insert Sizes:
if [[ $MODE == "PE" ]]; then
ls *_dedup.bam | \
  parallel "picard CollectInsertSizeMetrics I={} O={.}_InsertSizes.txt H={.}_InsertSizes.pdf QUIET=true VERBOSITY=ERROR 2> /dev/null"
  fi

####################################################################################################################################

## Library Complexity:
if [[ $MODE == "PE" ]]; then

  (>&2 paste -d " " <(echo '[INFO] LibComplexity started on') <(date))
  ls *_dup.bam | \
  parallel "preseq c_curve -bam -pe -s 5e+05 -o {.}_ccurve.txt {}"
  (>&2 paste -d " " <(echo '[INFO] LibComplexity ended on') <(date))
  fi
  
if [[ $MODE == "SE" ]]; then
  (>&2 paste -d " " <(echo '[INFO] LibComplexity started on') <(date))
  ls *_dup.bam | \
  parallel "preseq c_curve -bam -s 5e+05 -o {.}_ccurve.txt {}"
  (>&2 paste -d " " <(echo '[INFO] LibComplexity ended on') <(date))
  fi 

####################################################################################################################################

## Deprecated from v8 on, now use Genrich 8https://github.com/jsh58/Genrich/releases)
## Peaks for each sample:
## if [[ $GENOME == "mm10" ]]; then GFLAG="mm"; fi
## if [[ $GENOME == "hg38" ]]; then GFLAG="hs"; fi

## ls *_cutsites.bed.gz | awk -F "_cutsites.bed.gz" '{print $1}' | sort -u | \
##   parallel "$MACS callpeak -t {}_cutsites.bed.gz -n {} -g $GFLAG \
##                           --extsize 100 --shift -50 --nomodel --keep-dup=all -f BED --call-summits -q 0.01"

####################################################################################################################################

## Call peaks on individual replicates, excluding blacklisted regions and at 1% FDR:
(>&2 paste -d " " <(echo '[INFO] Peak Calling started on') <(date))

ls *_dedup.bam | awk -F "_dedup.bam" '{print $1}' | \
  parallel -j 8 "GenR {} 0.01 $BLACKLIST"
  
(>&2 paste -d " " <(echo '[INFO] Peak Calling ended on') <(date))

####################################################################################################################################

## FRiP:
(>&2 paste -d " " <(echo '[INFO] FRiP score calculation started on') <(date))

ls *cutsites.bed.gz | \
  awk -F "_cutsites.bed.gz" '{print $1}' | \
  parallel "FRiP {}_peaks_Genrich.narrowPeak {}_cutsites.bed.gz" > FRiP_scores.txt
  
(>&2 paste -d " " <(echo '[INFO] FRiP score calculation endedn on') <(date))

####################################################################################################################################

## Summary report:
multiqc -o multiqc_all ./ 
