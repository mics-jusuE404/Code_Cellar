#!/bin/bash
####################################################################################################################################
####################################################################################################################################
##
## Lowlevel script for automated processing of new mpra data.
## Naming conventions are ${BASENAME}_(R/D)NA_rep*.fastq.gz
## Steps:
#########
##    -- adapter trimming and alignment to hg38 to produce unfiltered, chr/quality-filtered & dup-filtered BAMs
##    -- calling peaks on pooled & deduplicated DNA samples
##    -- estimate library complexity and saturation with preseq, and make simply x/y plot
##    -- with the obtained peaks, make a count matrix for all dup- and dedup,
##    -- for this, first extend reads to fragment length (because 50bp single-end sequencing) and store as bed.gz
##    -- make CPM-normalized browser tracks for all BAM files
#########
## Written by Alexander Toenges (Nov 2018) a.toenges[(#aet#)]uni-muenster.de
## 
####################################################################################################################################
####################################################################################################################################
##
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=71
#SBATCH --partition=hims
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=mpra_lowlevel
#SBATCH --output=lowlevel.log
##
####################################################################################################################################
####################################################################################################################################

## Basic function for BAM quality check
function BamCheck {
  
  BASENAME=${1%_*}
  
  ## Run this to stop run in case of a suspiciously looking Bam:
  function ExitBam {
    echo '[ERROR]' $1 'looks suspicious or is empty -- exiting' 2>> ${BASENAME}.log 
  }; export -f ExitBam
  
  samtools quickcheck -q $1 && echo '' > /dev/null || ExitBam
  
  ## Also check if file is not empty:
  if [[ $(samtools view $1 | head -n 1 | wc -l) < 1 ]]; then
    ExitBam
    fi
  
}; export -f BamCheck  

## Adapter-trim and align data to hg38:
function Fq2Bam {

  BASENAME=$1
  
  echo '[Fq2Bam-START]' $BASENAME 'on:' >> ${BASENAME}.log && date >> ${BASENAME}.log && echo '' >> ${BASENAME}.log
  
  if [[ ! -e ${BASENAME}.fastq.gz ]]; then
    echo '[ERROR] Input file missing -- exiting' && exit 1; fi
  
  ## Adapters/Index:
  ADAPTER1="CTGTCTCTTATACACATCT"
  BWA_IDX=/scratch/tmp/a_toen03/Genomes/hg38/bwa_index_noALT_withDecoy/hg38_noALT_withDecoy.fa
  CHROMSIZES="/scratch/tmp/a_toen03/Genomes/hg38/hg38_chromSizes.txt"
  
  ## trim adapters, align and sort:
  cutadapt -j 4 -a $ADAPTER1 -m 36 --max-n 0.1 ${BASENAME}.fastq.gz 2>> ${BASENAME}.log | \
    bwa mem -v 2 \
      -R '@RG\tID:'${BASENAME}'_ID\tSM:'${BASENAME}'_SM\tPL:Illumina' -p -t 16 ${BWA_IDX} /dev/stdin 2>> ${BASENAME}.log | \
    samblaster --ignoreUnmated 2>> ${BASENAME}.log | \
    sambamba view -f bam -S -l 1 -t 4 -o /dev/stdout /dev/stdin 2>> ${BASENAME}.log | \
    sambamba sort -m 2G --tmpdir=./ -l 6 -t 16 -o ${BASENAME}_raw.bam /dev/stdin 2>> ${BASENAME}.log
    ##
    BamCheck ${BASENAME}_raw.bam
    
    
    ## kick out non-primary chromosomes and unmapped reads
    samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -vE 'chrM|_random|chrU|chrEBV|\*' 2>> ${BASENAME}.log | \
      xargs sambamba view -f bam -t 8 --num-filter=0/4 --filter='mapping_quality > 19' 2>> ${BASENAME}.log \
      -o ${BASENAME}_sortedDup.bam ${BASENAME}_raw.bam 2>> ${BASENAME}.log 
    ##
    BamCheck ${BASENAME}_sortedDup.bam
    
    ## kick out duplicated previously marked by samblaster
    sambamba view -f bam -t 8 --num-filter=/1028 -o ${BASENAME}_sortedDeDup.bam ${BASENAME}_sortedDup.bam 2>> ${BASENAME}.log
    ## 
    BamCheck ${BASENAME}_sortedDeDup.bam
    
    ## flagstat reports
    ls ${BASENAME}*.bam | parallel "sambamba flagstat -t 8 {} > {.}.flagstat 2>> ${BASENAME}.log"
    
  echo '[Fq2Bam-END]' $BASENAME 'on:' >> ${BASENAME}.log && date >> ${BASENAME}.log

}; export -f Fq2Bam

ls *.fastq.gz | awk -F ".fastq.gz" '{print $1}' | parallel -j 4 "Fq2Bam {}"

####################################################################################################################################
####################################################################################################################################

## Get the percentage of mitochondrial DNA in library using the samtools idxstats report:
function mtDNA {

  BASENAME=$1 
  
  BamCheck ${BASENAME}_raw.bam
  
  if [[ ! -e ${BASENAME}_raw.bam.bai ]]; then
    sambamba index ${BASENAME}_raw.bam
    fi
    
  mtReads=$(samtools idxstats ${BASENAME}_raw.bam | grep 'chrM' | cut -f 3)
  totalReads=$(samtools idxstats ${BASENAME}_raw.bam | awk '{SUM += $3} END {print SUM}')

 echo '[mtDNA Content]' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%' > ${BASENAME}_mtDNA.txt
 
}; export -f mtDNA

ls *.fastq.gz | awk -F ".fastq.gz" '{print $1}' | parallel "mtDNA {}"

####################################################################################################################################
####################################################################################################################################

## Combine all DNA/RNA files including duplicates into one file to run with preseq:
## will crash if one group has no replicates as sambamba expects at least 2 files for merging
function COMPLEXITY {
  
  BASENAME=$1

  ## Check involved BAMs:
  find ./ -maxdepth 1 -name "${BASENAME}*_rep*_sortedDup.bam" | \
    while read p 
      do; BamCheck $p; done < /dev/stdin
      
  ## Merge dup-BAMs for complexity check:
  find ./ -maxdepth 1 -name "${BASENAME}*_rep*_sortedDup.bam" | \
    xargs sambamba merge -t 16 ${BASENAME}_combined_sortedDup.bam
    
  BamCheck ${BASENAME}_combined_sortedDup.bam
  
  ## Run preseq c_curve to assess library complexity and saturation at given sequencing depth:
  find ./ -maxdepth 1 -name "*_sortedDup.bam" | \
    parallel "preseq c_curve -o {.}_ccurve -seed 1 -bam {}"
    
}; export -f COMPLEXITY

ls *.fastq.gz | awk -F "_rep" '{print $1}' | parallel -j 4 COMPLEXITY {}

####################################################################################################################################
####################################################################################################################################

function LosPeakos {
  
  BASENAME=$1
  
  ## Call peaks on the combined DNA sets per cell type with extsize corresponding to average fragment size,
  ## without q-value filtering, as this can be done manually afterwards on $9 of narrowPeak which is -log10(q)
  source activate py27
  
  while read p
    do; BamCheck $p
    done < <(ls ${BASENAME}*rep*sortedDeDup.bam)
    
  macs2 callpeak -t ${BASENAME}*rep*sortedDeDup.bam --nomodel --extsize 80 -n ${BASENAME} -g hs -f BAM
  
  ## Take highly-significant summits (q < 0.001), extend by 2 times the fragment size (160bp) and write as BED:
  bedtools slop -b 80 -g $CHROMSIZES -i ${BASENAME}_summits.bed | \
    sort -k1,1 -k2,2n > ${BASENAME}_referenceRegions.bed
    
}; export -f LosPeakos    

ls *.fastq.gz | awk -F "_rep" '{print $1}' | parallel LosPeakos {}

####################################################################################################################################
####################################################################################################################################

function CountMatrix {
  
  BASENAME=$1
  EXT=$2
  
  find ./ -maxdepth 1 -name "${BASENAME}*.bam" | grep -v 'combined' | \
    while read p
      do; BamCheck $p
      done < /dev/stdin
      
  ## First, write all BAM files as fragment-size (80bp) .bed.gz:
  function Bam2Bed {
    bedtools bamtobed -i $1 | \
      mawk -v LEN=${EXT} 'OFS="\t" {if ($6 == "+") {print $1, $2, $2+LEN} else if($6 == "-") print $1, $3-LEN, $3}' | \
      sort -k1,1 -k2,2n --parallel=4
  }; export -f Bam2Bed
 
  find ./ -maxdepth 1 -name "${BASENAME}*.bam" | grep -v 'combined' | \
    parallel "Bam2Bed {} | bgzip -@ 2 > {.}_ext${EXT}.bed.gz"
  
  ## Then, use bedtools intersect to write a count matrix with header:
  function IntersectX {

    cat <(ls ${BASENAME}*_ext.bed.gz | xargs echo 'chr' 'start' 'end' | tr " " "\t") \
        <(bedtools intersect -sorted -c -a ${BASENAME}_referenceRegions.bed -b ${BASENAME}_ext${EXT}.bed.gz | cut -f1-3,7-)
        
  }; export -f IntersectX
    
  IntersectX ${BASENAME}
}; export -f CountMatrix    

## CountMatrix Basename FragmentSize
ls *.fastq.gz | awk -F "_rep" '{print $1}' | parallel "CountMatrix {} 80"

####################################################################################################################################
####################################################################################################################################

## Bigwigs:
ls *_.bam | grep -vE 'raw' | parallel -j 8 "bamCoverage -e 80 --normalizeUsing CPM -bs 1 --bam {} -o {.}_e80_CPM.bigwig -p 9" 
