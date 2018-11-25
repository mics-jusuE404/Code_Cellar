#!/bin/bash
########
########
## Lowlevel processing for mpra data.
## Naming conventions for fastqs are ${BASENAME}_(R/D)NA_rep*.fastq.gz
#########
#########
## Steps:
##    -- adapter trimming and alignment to hg38 to produce unfiltered, chr/quality-filtered & dup-filtered BAMs
##    -- calling peaks on pooled & deduplicated DNA samples
##    -- estimate library complexity and saturation with preseq, and make simply x/y plot
##    -- with the obtained peaks, make a count matrix for all dup- and dedup,
##    -- for this, first extend reads to fragment length (because 50bp single-end sequencing) and then count reads over reference
##    -- make CPM-normalized browser tracks for all BAM files
#########
#########
##    -- Written by Alexander Toenges (Nov 2018) a.toenges[(#aet#)]uni-muenster.de 
#########
#########
####################################################################################################################################
####################################################################################################################################
###
##
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=71
#SBATCH --partition=hims
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=mpra_lowlevel
#SBATCH --output=lowlevel.log
#
##
###
####################################################################################################################################
####################################################################################################################################

CHROMSIZES="/scratch/tmp/a_toen03/Genomes/hg38/hg38_chromSizes.txt"
ADAPTER1="CTGTCTCTTATACACATCT"
BWA_IDX="/scratch/tmp/a_toen03/Genomes/hg38/bwa_index_noALT_withDecoy/hg38_noALT_withDecoy.fa"

####################################################################################################################################
####################################################################################################################################

## Quality check:
function ExitBam {

  (>&2 echo '[ERROR]' $1 'looks suspicious or is empty -- exiting') && exit 1
  
}; export -f ExitBam

####################################################################################################################################
####################################################################################################################################

## Basic function for BAM quality check
function BamCheck {
  
  BASENAME=${1%_*}
  samtools quickcheck -q $1 && echo '' >/dev/null || ExitBam $1
  
  ## Also check if file is not empty:
  if [[ $(samtools view $1 | head -n 1 | wc -l) < 1 ]]; then
    ExitBam $BASENAME
    fi
  
}; export -f BamCheck  

####################################################################################################################################
####################################################################################################################################

## Adapter-trim and align data to hg38:
function Fq2Bam {
  
  BASENAME=$1
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'started on') <(date))
  
  if [[ ! -e ${BASENAME}.fastq.gz ]]; then
    (>&2 echo '[ERROR] Input file missing -- exiting') && exit 1; fi
  
  ## trim adapters, align and sort:
  cutadapt -j 4 -a $ADAPTER1 -m 36 --max-n 0.1 ${BASENAME}.fastq.gz | \
    bwa mem -v 2 -R '@RG\tID:'${BASENAME}'_ID\tSM:'${BASENAME}'_SM\tPL:Illumina' -p -t 16 ${BWA_IDX} /dev/stdin | \
    samblaster --ignoreUnmated | \
    sambamba view -f bam -S -l 1 -t 4 -o /dev/stdout /dev/stdin | \
    sambamba sort -m 2G --tmpdir=./ -l 6 -t 16 -o ${BASENAME}_raw.bam /dev/stdin
    ##
    BamCheck ${BASENAME}_raw.bam
    
    ## kick out non-primary chromosomes and unmapped reads
    samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -vE 'chrM|_random|chrU|chrEBV|\*' | \
      xargs sambamba view -f bam -t 8 --num-filter=0/4 --filter='mapping_quality > 19' \
      -o ${BASENAME}_sortedDup.bam ${BASENAME}_raw.bam
    ##
    BamCheck ${BASENAME}_sortedDup.bam
    
    ## kick out duplicated previously marked by samblaster
    sambamba view -f bam -t 8 --num-filter=/1028 -o ${BASENAME}_sortedDeDup.bam ${BASENAME}_sortedDup.bam
    ## 
    BamCheck ${BASENAME}_sortedDeDup.bam
    
    ## flagstat reports
    ls ${BASENAME}*.bam | parallel "sambamba flagstat -t 8 {} > {.}.flagstat"
    
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'ended on') <(date))

}; export -f Fq2Bam

ls *.fastq.gz | awk -F ".fastq.gz" '{print $1}' | parallel -j 4 "Fq2Bam {} 2>> {}.log"

####################################################################################################################################
####################################################################################################################################

## Get the percentage of mitochondrial DNA in library using the samtools idxstats report:
function mtDNA {
  
  (>&2 paste -d " " <(echo '[INFO]' 'mtDNA for' $1 'started on') <(date))
  
  BASENAME=$1 
  
  BamCheck ${BASENAME}_raw.bam
  
  if [[ ! -e ${BASENAME}_raw.bam.bai ]]; then
    sambamba index ${BASENAME}_raw.bam
    fi
    
  mtReads=$(samtools idxstats ${BASENAME}_raw.bam | grep 'chrM' | cut -f 3)
  totalReads=$(samtools idxstats ${BASENAME}_raw.bam | awk '{SUM += $3} END {print SUM}')

 echo '[mtDNA Content]' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%' > ${BASENAME}_mtDNA.txt
 
 (>&2 paste -d " " <(echo '[INFO]' 'mtDNA for' $1 'ended on') <(date))
 
}; export -f mtDNA

ls *.fastq.gz | awk -F ".fastq.gz" '{print $1}' | parallel "mtDNA {} 2>> {}.log"

####################################################################################################################################
####################################################################################################################################

## Combine all _sortedDup.bam files into one for each condition for use with preseq c_curve:
function MergeBam {
  
  (>&2 paste -d " " <(echo '[INFO]' 'MergeBam for' $1 'started on') <(date))
  
  BASENAME=$1

  find ./ -maxdepth 1 -name "${BASENAME}*_rep*_sortedDup.bam" | \
    while read p 
      do BamCheck $p; done < /dev/stdin
      
  find ./ -maxdepth 1 -name "${BASENAME}*_rep*_sortedDup.bam" | \
    xargs sambamba merge -t 16 ${BASENAME}_merged_sortedDup.bam
  
  (>&2 paste -d " " <(echo '[INFO]' 'MergeBam for' $1 'ended on') <(date))
  
}; export -f MergeBam

ls *_sortedDup.bam | awk -F "_rep" '{print $1}' | sort -k1,1 -u | grep -v 'RNA' | \
  parallel -j 4 "MergeBam {} 2>> {.}_MergeBam.log"

## Combine all DNA/RNA files including duplicates into one file to run with preseq:
## will crash if one group has no replicates as sambamba expects at least 2 files for merging

####################################################################################################################################
####################################################################################################################################

function Complexity {
  
  (>&2 paste -d " " <(echo '[INFO]' 'Complexity for' $1 'started on') <(date))

  BAM=$1
  BamCheck $BAM
  
  ## Create library complexity curve:
  preseq c_curve -s 5e+05 -o ${BAM%.bam}_ccurve.txt -seed 1 -bam $BAM
  
  (>&2 paste -d " " <(echo '[INFO]' 'COMPLEXITY for' $1 'ended on') <(date))
    
}; export -f Complexity

ls *DNA*.bam | grep -vE 'DeDup|raw' | parallel "Complexity {} {.}_complexity.log"

####################################################################################################################################
####################################################################################################################################

function PeakCall {
  
  BASENAME=$1
  
  (>&2 paste -d " " <(echo '[INFO]' 'PeakCall for' $1 'started on') <(date))
  
  ## Call peaks on the combined DNA sets per cell type with extsize corresponding to average fragment size,
  ## without q-value filtering, as this can be done manually afterwards on $9 of narrowPeak which is -log10(q)
  
  ls ${BASENAME}*rep*sortedDeDup.bam | while read p
    do BamCheck $p
    done < /dev/stdin
  
  source activate py27
  macs2 callpeak --nomodel --extsize 80 -g hs -f BAM -n ${BASENAME} -t ${BASENAME}*rep*sortedDeDup.bam
  source deactivate
  
  ## Take highly-significant summits (q < 0.001), extend by 2 times the fragment size (160bp) and write as BED:
  bedtools slop -b 80 -g "$CHROMSIZES" -i ${BASENAME}_summits.bed | \
    sort -k1,1 -k2,2n > ${BASENAME}_referenceRegions.bed
  
  awk 'OFS="\t" {print $1":"$2+1"-"$3, $1, $2+1, $3, "+"}' ${BASENAME}_referenceRegions.bed > ${BASENAME}_referenceRegions.saf  
  
  (>&2 paste -d " " <(echo '[INFO]' 'PeakCall for' $1 'ended on') <(date))

}; export -f PeakCall    

ls *sortedDeDup.bam | awk -F "_rep" '{print $1}' | sort -k1,1 -u | grep -v 'RNA' | parallel "LosPeakos {} 2>> {}_macs2.log"

####################################################################################################################################
####################################################################################################################################

## Extend all reads to fragment length:
function BamExtend {
  
  (>&2 paste -d " " <(echo '[INFO]' 'BamExtend for' $1 'started on') <(date))
  
  BAM=$1
  EXT=$2
  
  find ./ -maxdepth 1 -name "${BASENAME}*.bam" | grep -v '_raw' | \
    while read p
      do BamCheck $p
      done < /dev/stdin
      
  bedtools bamtobed -i $BAM | \
    mawk -v ext=${EXT} 'OFS="\t" {if ($6 == "+") {print $1,$2,$2+ext,$4,$5,$6} else if($6 == "-") print $1,$3-ext,$3,$4,$5,$6}' | \
    bedtools bedtobam -i - -g $CHROMSIZES > ${BAM%.bam}_ext${EXT}.bam
    
  (>&2 paste -d " " <(echo '[INFO]' 'BamExtend for' $1 'ended on') <(date))
  
}; export -f BamExtend
  
ls *.bam | grep -vE 'merged|_raw' | \
  parallel "BamExtend {} 80 2>> {.}_BamExtend.log"

####################################################################################################################################
####################################################################################################################################

function CountMatrix {
  
  BASENAME=$1
    
  ## Then, use bedtools intersect to write a count matrix with header:
  function Fcount {

    featureCounts -T 8 -a ${BASENAME}_referenceRegions.saf -o ${BASENAME}.counts \
      ${BASENAME}*rep*sortedDup.bam ${BASENAME}*rep*sortedDeDup.bam
        
  }; export -f Fcount
    
  Fcount ${BASENAME}
}; export -f CountMatrix    

## CountMatrix Basename FragmentSize
ls *.fastq.gz | awk -F "_rep" '{print $1}' | sort -k1,1 -u | \
  'parallel "CountMatrix {} 2>> {}_countmatrix.log"

####################################################################################################################################
####################################################################################################################################

## Bigwigs:
ls *_.bam | grep -vE 'raw' | parallel -j 8 "bamCoverage -e 80 --normalizeUsing CPM -bs 1 --bam {} -o {.}_e80_CPM.bigwig -p 9" 
