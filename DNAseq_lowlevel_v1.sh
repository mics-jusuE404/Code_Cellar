#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --partition=normal
#SBATCH --mem=80G
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=DNAseq_lowlevel

######################################################################################################################################
######################################################################################################################################

## Script for lowlevel processing of DNA-seq data: Alignment, filtering, quality control, browser tracks, peak calling, 
## FRiPs, library complexity, multiqc summary etc...

######################################################################################################################################
######################################################################################################################################

## Customizable parameters:
## => MODE can be UBAMSE/UBAMPE/FQSE/FQPE
## => JOBS can be set to 8 if using the SMP nodes or 4 if on the normal Skylake nodes
## => ASSAY can be ATACseq or ChIPseq

ASSAY="ATACseq"
MODE="UBAMPE"
GENOME="mm10"
JOBS=4

######################################################################################################################################
######################################################################################################################################

RSCRIPT="$HOME/anaconda3_things/anaconda3/envs/R_env/bin/Rscript"
MACS="$HOME/anaconda3_things/anaconda3/envs/py27_env/bin/macs2"

if [[ ${GENOME} == "hg38" ]]; then
  IDX="/scratch/tmp/a_toen03/Genomes/hg38/bowtie2_index_noALT_withDecoy/hg38_noALT_withDecoy.fa"
  BLACKLIST="/scratch/tmp/a_toen03/Genomes/hg38/Blacklists/hg38_consensusBL.bed"
  fi
  
if [[ ${GENOME} == "mm10" ]]; then
  IDX="/scratch/tmp/a_toen03/Genomes/mm10/bowtie2_idx/mm10"
  BLACKLIST="/scratch/tmp/a_toen03/Genomes/mm10/mm10_consensusBL.bed"
  fi  

######################################################################################################################################
######################################################################################################################################

## => Section to check parameters and tools being in PATH:

## Check if MODE and GENOME are correctly set is correct:
if [[ ! $(echo $MODE | grep -E 'UBAMSE|UBAMPE|FQSE|FQPE') ]]; then 
  '[ERROR] MODE parameter is neither of (UBAMSE|UBAMPE|FQSE|FQPE) -- exiting' && exit 1
  fi
  
if [[ ! $(echo $GENOME | grep -E 'hg38|mm10') ]]; then 
  '[ERROR] GENOME is neither of the supported hg38 and mm10 -- exiting' && exit 1
  fi

if [[ ! $(echo $MODE | grep -E 'ATACseq|ChIPseq') ]]; then 
  '[ERROR] ASSAY is neither of the supported ATACseq or ChIPseq -- exiting' && exit 1
  fi

## Check if required tools are in PATH:
if [[ -e missing_tools.txt ]]; then rm missing_tools.txt; fi

## Function that checks if required tools are callable from PATH:
function PathCheck {
  if [[ $(command -v $1 | wc -l) == 0 ]]; then 
    echo ${1} >> missing_tools.txt
    fi
}; export -f PathCheck

## All tools:
TOOLS=(cat seqtk cutadapt bwa samtools samblaster sambamba tee xargs bedtools mawk bgzip tabix \
       sort paste featureCounts bc bamCoverage parallel fastqc picard preseq multiqc $RSCRIPT \
       bigWigToBedGraph $MACS bg2bw pigz featureCounts)

## Loop through tools:
for i in $(echo ${TOOLS[*]}); do
  PathCheck $i
  done
  
## If tools are missing write them to <missing_tools.txt>
if [[ -e missing_tools.txt ]] && [[ $(cat missing_tools.txt | wc -l | xargs) > 0 ]]; then
  echo '[ERROR] Tools missing in PATH -- see missing_tools.txt' && exit 1
  fi

######################################################################################################################################
######################################################################################################################################

## Rscript to calculate size factors using DESeq2 which is printed to disk:
cat <<EOF > sizeFactors.R
## Script to calculate DESeq2 size factors for ATAC-seq data based on a count matrix
## over the entire genome with 500bp windows, taking the top 100k windows based on rowMeans.
## File comes from stdin via the main bash script ATACseq_lowlevel(...).sh

packageS <- c("DESeq2", "data.table")
  if (length(grep("FALSE", (packageS %in% rownames(installed.packages())))) > 0){
    stop("Package(s): ", packageS[which( packageS %in% rownames(installed.packages()) == "FALSE")], " are not installed!")
}

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
counts <- data.frame(counts, rowMean = rowMeans(counts[,7:ncol(counts)]), check.names=F)

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
######################################################################################################################################

## Exit function if BAM file looks corrupted or is missing:
function ExitBam {

  (>&2 echo '[ERROR]' "${1}" 'looks suspicious or is empty -- exiting') && exit 1
  
}; export -f ExitBam

######################################################################################################################################
######################################################################################################################################

## Check if BAM is corrupted or empty:
function BamCheck {
  
  BASENAME="${1%_*}"
  samtools quickcheck -q $1 && echo '' >/dev/null || ExitBam "${1}"
    
  ## Also check if file is not empty:
  if [[ $(samtools view "${1}" | head -n 1 | wc -l) < 1 ]]; then
    ExitBam $BASENAME
    fi
  
}; export -f BamCheck  

######################################################################################################################################
######################################################################################################################################

## Function to get the % of reads mapped to chrM:
function mtDNA {

  BamCheck $1
    
  if [[ ! -e "${1}".bai ]];
    then
    sambamba index -t 8 ${1}
    fi
      
  mtReads=$(samtools idxstats $1 | grep 'chrM' | cut -f 3)
  totalReads=$(samtools idxstats $1 | awk '{SUM += $3} END {print SUM}')

  echo '[mtDNA Content]' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%' > "${1%.bam}"_mtDNA.txt
  
}; export -f mtDNA

######################################################################################################################################
######################################################################################################################################

function Fq2Bam {

  BASENAME="${1}"
  TYPE="${2}"
  IDX="${3}"
  ASSAY="${4}"
	
  #############################################################################################################
  ## 4 versions upstream of bowtie2 depending on paired/single and fq or ubam:
  if [[ ${ASSAY} == "ATACseq" ]]; then
    ADAPTER="CTGTCTCTTATACACATCT"; fi
  
  if [[ ${ASSAY} == "ChIPseq" ]]; then
    ADAPTER="AGATCGGAAGAGC"; fi  
    
  CUT_FQSE="cutadapt -j 1 -a "${ADAPTER}" -m 18 --max-n 0.1 "${BASENAME}".fastq.gz"
	
  CUT_FQPE="seqtk mergepe "${BASENAME}"_1.fastq.gz "${BASENAME}"_2.fastq.gz \
  	    | cutadapt -j 1 -a "${ADAPTER}" -A "${ADAPTER}" --interleaved -m 18 --max-n 0.1 -"
  	        
  CUT_UBSE="samtools fastq -n -@ 2 "${BASENAME}"_ubam.bam \
            | cutadapt -j 1 -a "${ADAPTER}" -m 18 --max-n 0.1 -"	 
            
  CUT_UBPE="samtools fastq -n -@ 2 "${BASENAME}"_ubam.bam \
            | cutadapt -j 1 -a "${ADAPTER}" -A "${ADAPTER}" --interleaved -m 18 --max-n 0.1 -"
  
  #############################################################################################################          
  
  ## two bowtie2 modes:
  BOWTIE2PE="bowtie2 --very-sensitive --threads 16 -X 2000 --rg-id "${BASENAME}" -x $IDX --interleaved - \
             | samtools fixmate -m -@ 2 -O SAM - -"

  BOWTIE2SE="bowtie2 --very-sensitive --threads 16 -x "${IDX}" -U -"             
	
  #############################################################################################################
  ## Check mode and set concat command accordingly:
  
  ## 1) uBAM-SE:
  if [[ $(echo ${TYPE} | grep 'UBAMSE') ]]; then
    MESSAGE1="UBAM-SE"
    COMMAND1="${CUT_UBSE} \
              | ${BOWTIE2SE}"; fi
  ## 2) uBAM-PE:
  if [[ $(echo ${TYPE} | grep 'UBAMPE') ]]; then    
    MESSAGE1="UBAM-PE"
    COMMAND1="${CUT_UBPE} \
              | ${BOWTIE2PE}"; fi
  ## 3) FQ-SE:              
  if [[ $(echo ${TYPE} | grep 'FQSE') ]]; then    
    MESSAGE1="FASTQ-SE"
    COMMAND1="${CUT_FQSE} \
              | ${BOWTIE2SE}"; fi
  ## 4) uBAM-PE:
  if [[ $(echo ${TYPE} | grep 'FQPE') ]]; then 
    MESSAGE1="FASTQ-PE"
    COMMAND1="${CUT_FQPE} \
              | ${BOWTIE2PE}"; fi              
  
  ## Summary message:    
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' "${1}" 'in' ${MESSAGE1} 'started on') <(date))

  ## Run the lowlevel alignment pipeline:
  eval "${COMMAND1}" \
  | samblaster --ignoreUnmated \
  | sambamba view -f bam -S -l 5 -t 2 -o /dev/stdout /dev/stdin \
  | tee >(sambamba flagstat -t 2 /dev/stdin > "${BASENAME}"_raw.flagstat) \
  | sambamba sort -m 5G --tmpdir=./ -l 5 -t 16 -o "${BASENAME}"_raw.bam /dev/stdin  
  
  ## Check if alignment went ok based on BAM file integrity:
  BamCheck "${BASENAME}"_raw.bam
  
  ## Get mtDNA percentage:
  mtDNA "${BASENAME}"_raw.bam
  
  ## Start filtering the BAM file:
  if [[ ! -e tmp_chromSizes.txt ]]; then
    samtools idxstats "${BASENAME}"_raw.bam > tmp_chromSizes.txt
    fi
  
  ## Decide appropriate flag to keep only mapped reads depending on paired/single-end sequencing:
  if [[ $(echo ${TYPE} | grep 'PE') ]]; then
    NUMFILTER=1
    fi
  if [[ $(echo ${TYPE} | grep 'SE') ]]; then
    NUMFILTER=0
    fi
	
	## Now define filtering that happens for all DNA-seq:
	GENERALFILTER="samtools idxstats "${BASENAME}"_raw.bam \
                 | cut -f 1 \
                 | grep -vE 'chrM|_random|chrU|chrEBV|\*' \
                 | xargs sambamba view -l 5 -f bam -t 8 --num-filter="${NUMFILTER}"/2308 --filter='mapping_quality > 19' \
                   -o /dev/stdout "${BASENAME}"_raw.bam \
                 | tee "${BASENAME}"_dup.bam \
                 | tee >(samtools index - "${BASENAME}"_dup.bam.bai)"
	
	## Specific for ATAC-seq (cutting site extraction):
	ATACspecific="sambamba view -l 5 -f bam -t 8 --num-filter=/1024 -o /dev/stdout /dev/stdin \
                | tee >(tee "${BASENAME}"_dedup.bam | samtools index - "${BASENAME}"_dedup.bam.bai) \
                | bedtools bamtobed -i - \
                | mawk 'OFS="\t" {if ($6 == "+") print $1, $2+4, $2+5, ".", ".", $6} {if ($6 == "-") print $1, $3-5, $3-4, ".", ".", $6}' \
                | sort -k1,1 -k2,2n -k3,3n -k6,6 -S10G --parallel=10 \
                | tee >(bgzip -@ 6 > "${BASENAME}"_cutsites.bed.gz) \
                | bedtools genomecov -bg -i - -g tmp_chromSizes.txt \
                | bg2bw -i /dev/stdin -c tmp_chromSizes.txt -o "${BASENAME}"_cutsites_noScale.bigwig"
  
  ## and non-ATACseq:              
  NONATAC="sambamba view -l 5 -f bam -t 8 --num-filter=/1024 -o ${BASENAME}_dedup.bam /dev/stdin"                
  
  ## Cat commands together:            	
  if [[ ${ASSAY} == "ChIPseq" ]]; then
    eval "${GENERALFILTER}" \
    | eval "${NONATAC}"
    fi
  
  if [[ ${ASSAY} == "ATACseq" ]]; then
    eval "${GENERALFILTER}" \
    | eval "${ATACspecific}"
    fi
  
  ## Check if BAMs look ok:
  BamCheck "${BASENAME}"_dup.bam
  BamCheck "${BASENAME}"_dedup.bam
  
  ## Flagstat that I couldn't squeeze into the tee part above:
  ls "${BASENAME}"*dup.bam | parallel "sambamba flagstat -t 8 {} > {.}.flagstat"
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' "${1}" 'in' "${MESSAGE1}" 'ended on') <(date))

}

export -f Fq2Bam

######################################################################################################################################

## Function that calculated DESeq2 size factors based on the 100k 500bp windows across the genome with the highest counts:
function SizeFactor {
  
  (>&2 paste -d " " <(echo '[INFO] SizeFactor started on') <(date))

  ## 500bp windows over the genome:
  bedtools makewindows -w 500 -g tmp_chromSizes.txt \
    | mawk 'OFS="\t" {print $1"_"$2"_"$3, $1, $2+1, $3, "+"}' > genome_windows.saf
     
  ## count matrix
  featureCounts --read2pos 5 -a genome_windows.saf -F SAF -T 8 -o genome_windows_counts.txt *_dedup.bam
  
  ## feed counts into the DESeq2 script to get normalized counts: 
  cat genome_windows_counts.txt \
    | "${RSCRIPT}" sizeFactors.R
  
  (>&2 paste -d " " <(echo '[INFO] SizeFactor ended on') <(date))
  
}; export -f SizeFactor
   
####################################################################################################################################

## Function to use the above size factors for bigwig tracks:
function Bigwig_SF {

  FILE="${1}"
    
  (>&2 paste -d " " <(echo '[INFO]' 'Browser Tracks for' "${1}" 'started on') <(date))
  
  ## bamcoverage multiplies with the size factor but deseq2 divides so invert the deseq factor:
  if [[ $(cat sizeFactors.txt | wc -l | xargs) > 1 ]]; then
    FACTOR=$(grep ${FILE%_cutsites*}_dedup.bam sizeFactors.txt \
      | cut -f2 \
      | bc <<< "scale=6;$(grep ${FILE%_cutsites*}_dedup.bam sizeFactors.txt | cut -f2)^-1")
    fi
    
  if [[ $(cat sizeFactors.txt | wc -l | xargs) == 1 ]]; then
    FACTOR=1
    fi
  
  ## Bigwig with 80bp fragments to smooth the signal:
  bigWigToBedGraph "${FILE}" /dev/stdout \
  | bedtools slop -b 40 -g tmp_chromSizes.txt -i - \
  | bedtools genomecov -g tmp_chromSizes.txt -i - -bg -scale "${FACTOR}" \
  | bg2bw -i /dev/stdin -c tmp_chromSizes.txt -o "${FILE%_noScale.bigwig}"_geomean.bigwig
  
  (>&2 paste -d " " <(echo '[INFO]' 'Browser Tracks for' "${1}" 'ended on') <(date))
  
}; export -f Bigwig_SF

## ## Function to make boring count-per-million bigwigs 
function Bigwig_CPM {

  FILE="${1}"
  MODE="${2}"
  
  (>&2 paste -d " " <(echo '[INFO]' 'Browser Tracks for' "${1}" 'started on') <(date))
  
  if [[ $(echo $MODE | grep 'PE') ]]; then
    Ext="-e"; fi
  
  if [[ $(echo $MODE | grep 'SE') ]]; then
    Ext="-e 150"; fi  
  
  bamCoverage --bam $FILE -o ${FILE%_dedup.bam}_CPM.bigwig ${Ext} -p 16 --normalizeUsing CPM
    
  (>&2 paste -d " " <(echo '[INFO]' 'Browser Tracks for' "${1}" 'ended on') <(date))  
  
}; export -f Bigwig_CPM

####################################################################################################################################

## FRiPs based on the cutting site bigwigs for ATAC-seq:
function FRiP_ATAC {
  
  PEAKS="${1}"
  CUTSITES="${2}"
  
  if [[ -e "${1}" ]] && [[ -e "${2}" ]]; then
    TOTAL=$(pigz -c -d -p 4 "${CUTSITES}" | wc -l)
    READ=$(bedtools intersect \
           -a "${CUTSITES}" \
           -b "${PEAKS}" -u \
           -wa -sorted | wc -l)
     paste <(echo "${2%_cutsites.bed.gz}") <(bc <<< "scale=6;$READ/$TOTAL") 
     fi

}; export -f FRiP_ATAC

## FRiP for nopn-ATAC based on the BAM files using featureCounts:
function FRiP_ChIP {
  
  PEAKS=$1
  BAM=$2
  MODE=$3
  
  if [[ -e $1 ]] && [[ -e $2 ]]; then
  
    if [[ ${MODE} == "PE" ]]; then
      INS="-p"
      fi
    if [[ ${MODE} == "SE" ]]; then
      INS=""
      fi 
      
    ## Peak to SAF:
    awk 'OFS="\t" {print $1$2$3, $1, $2+1, $3, "+"}' ${PEAKS} > ${PEAKS}.saf
    
    ## Count with featureCounts:
    featureCounts "${INS}" -T 16 -a ${PEAKS}.saf -F SAF -o ${PEAKS}_fc.txt $BAM
    
    ASSIGNED=$(grep -w 'Assigned' ${PEAKS}_fc.txt.summary | cut -f2)
    UNASSIGNED=$(grep -w 'Unassigned_NoFeatures' ${PEAKS}_fc.txt.summary | cut -f2)
    
    rm ${PEAKS}_fc.txt ${PEAKS}_fc.txt.summary
    
    paste <(echo "${PEAKS%_noControl*}") <(bc <<< "scale=6;${ASSIGNED}/(${ASSIGNED}+${UNASSIGNED})")
    fi

  }; export -f FRiP_ChIP

####################################################################################################################################
####################################################################################################################################
##
## Call above functions:
##
####################################################################################################################################
####################################################################################################################################

## fastqc:
if [[ $(echo $MODE | grep 'UBAM') ]]; then
  ls *_ubam.bam | parallel -j 70 "fastqc -t 1 {}"; fi
  
if [[ $(echo $MODE | grep 'FQ') ]]; then
  ls *.fastq.gz | parallel -j 70 "fastqc -t 1 {}"; fi  
  
####################################################################################################################################

## Alignment/filtering function:
if [[ $(echo $MODE | grep 'UBAM') ]]; then
  ls *_ubam.bam | awk -F "_ubam.bam" '{print $1}' | sort -u | parallel -j $JOBS "Fq2Bam {} ${MODE} ${IDX} ${ASSAY} 2>> {}.log"
  fi

if [[ $(echo $MODE | grep 'FQPE') ]]; then
  ls *_1.fastq.gz | awk -F "_1.fastq.gz" '{print $1}' | sort -u | parallel -j $JOBS "Fq2Bam {} ${MODE} ${IDX} ${ASSAY} 2>> {}.log"
  fi

if [[ $(echo $MODE | grep 'FQSE') ]]; then
  ls *.fastq.gz | awk -F ".fastq.gz" '{print $1}' | sort -u | parallel -j $JOBS "Fq2Bam {} ${MODE} ${IDX} ${ASSAY} 2>> {}.log"
  fi

####################################################################################################################################

## Estimate size factors using DESeq2 but only for ATAC-seq
## (for ChIP-seq pretty meaningless given that input samples would skew the analysis):
if [[ $ASSAY == "ATACseq" ]]; then

  if [[ $(ls *_raw.bam | wc -l) > 1 ]]; then
  
    SizeFactor 2> sizeFactors.log
    else  
  
    echo 'NULL' > sizeFactors.txt
    echo 'Only one sample' > sizeFactors.log
  
    fi
  fi  
  
####################################################################################################################################

## Get browser tracks, scaled by the size factor from DESeq2 (but only if more than one sample present,
## else would be meaningless:

if [[ ${ASSAY} == "ATACseq" ]]; then

  if [[ $(ls *_cutsites_noScale.bigwig | wc -l) > 1 ]]; then
  
    ls *_cutsites_noScale.bigwig \
      | awk -F "_cutsites_noScale.bigwig" '{print $1}' \
      | parallel -j 16 "Bigwig {}_cutsites_noScale.bigwig 2>> {}.log"
    else
      (>&2 echo '[INFO] Only one sample present, skipping normalization of that one bigwig file')
    fi
  fi

## ChIP-seq tracks just simple CPM:  
if [[ ${ASSAY} == "ChIPseq" ]]; then

  ls *_dedup.bam | parallel -j $JOBS "Bigwig_CPM {} ${MODE}"
  
  fi
  
####################################################################################################################################

## Insert Sizes given paired-end data:
if [[ $(echo $MODE | grep 'PE') ]]; then

  (>&2 paste -d " " <(echo '[INFO] CollectInsertSizes started on') <(date))
  
  ls *_dedup.bam \
    | parallel "picard CollectInsertSizeMetrics I={} O={.}_InsertSizes.txt H={.}_InsertSizes.pdf \
                QUIET=true VERBOSITY=ERROR VALIDATION_STRINGENCY=LENIENT 2> /dev/null"
	    
  (>&2 paste -d " " <(echo '[INFO] CollectInsertSizes ended on') <(date))
  
  fi

####################################################################################################################################

## Library Complexity:
if [[ $(echo $MODE | grep 'PE') ]]; then
  (>&2 paste -d " " <(echo '[INFO] LibComplexity started on') <(date))
  ls *_dup.bam | \
  parallel "preseq c_curve -bam -pe -s 5e+05 -o {.}_ccurve.txt {}"
  (>&2 paste -d " " <(echo '[INFO] LibComplexity ended on') <(date))
  fi
  
if [[ $(echo $MODE | grep 'SE') ]]; then
  (>&2 paste -d " " <(echo '[INFO] LibComplexity started on') <(date))
  ls *_dup.bam | \
  parallel "preseq c_curve -bam -s 5e+05 -o {.}_ccurve.txt {}"
  (>&2 paste -d " " <(echo '[INFO] LibComplexity ended on') <(date))
  fi 

####################################################################################################################################

## Peak calling with macs:
if [[ $GENOME == "mm10" ]]; then GFLAG="mm"; fi
if [[ $GENOME == "hg38" ]]; then GFLAG="hs"; fi

## change the name of the peaks in the filtered file from unfiltered to filtered:		
function GSUB {
	
  awk 'OFS="\t" {gsub("unfiltered","filtered");print}' "${1}" > "${1}".tmp
  mv "${1}".tmp "${1}"
  
}; export -f GSUB

## => 
## Call Peaks for ATAC-seq at 1% FDR:
if [[ ${ASSAY} == "ATACseq" ]]; then

  ls *_cutsites.bed.gz \
    | awk -F "_cutsites.bed.gz" '{print $1}' \
    | sort -u \
    | parallel "$MACS callpeak -t {}_cutsites.bed.gz -n {}_FDR1perc_unfiltered -g $GFLAG \
                --extsize 100 --shift -50 --nomodel --keep-dup=all -f BED --call-summits -q 0.01 \
  		          --min-length 150"
  		          
  ## Filter against blacklist:
  ls *_FDR1perc_unfiltered_peaks.narrowPeak \
    | awk -F "_FDR1perc_unfiltered_peaks.narrowPeak" '{print $1}' \
    | parallel "bedtools intersect -v -a {}_FDR1perc_unfiltered_peaks.narrowPeak \
                                   -b ${BLACKLIST} > {}_FDR1perc_filtered_peaks.narrowPeak"
  				  
  ls *_FDR1perc_unfiltered_summits.bed \
    | awk -F "_FDR1perc_unfiltered_summits.bed" '{print $1}' \
    | parallel "bedtools intersect -v -a {}_FDR1perc_unfiltered_summits.bed \
                                   -b ${BLACKLIST} > {}_FDR1perc_filtered_summits.bed"
  
  find -maxdepth 1 -name "*_filt*" -printf '%f\n' \
  | parallel "GSUB {}"
  
  fi
  
## ChIP-seq at 5% FDR
if [[ ${ASSAY} == "ChIPseq" ]]; then

  ls *_dedup.bam \
    | awk -F "_dedup.bam" '{print $1}' \
    | sort -u \
    | parallel "$MACS callpeak -t {}_dedup.bam -n {}_noControl_unfiltered -g $GFLAG \
                --keep-dup=all"
  		          
  ## Filter against blacklist:
  ls *_FDR5perc_unfiltered_peaks.narrowPeak \
    | awk -F "_noControl_unfiltered_peaks.narrowPeak" '{print $1}' \
    | parallel "bedtools intersect -v -a {}_noControl_unfiltered_peaks.narrowPeak \
                                   -b ${BLACKLIST} > {}_noControl_filtered_peaks.narrowPeak"
  				  
  find -maxdepth 1 -name "*_filt*" -printf '%f\n' \
  | parallel "GSUB {}"
  
  fi
  
####################################################################################################################################

## Get Fractions of Reads per Peak (FRiPs):
(>&2 paste -d " " <(echo '[INFO] FRiP score calculation started on') <(date))

if [[ ${ASSAY} == "ATACseq" ]]; then
  
  ls *cutsites.bed.gz \
    | awk -F "_cutsites.bed.gz" '{print $1}' \
    | parallel "FRiP_ATAC {}_FDR1perc_filtered_peaks.narrowPeak {}_cutsites.bed.gz" > FRiP_scores.txt
    
  fi  

## FRiPs for ChIP-seq from BAM files  
if [[ ${ASSAY} == "ChIPseq" ]]; then

  ls *_noControl_peaks.narrowPeak \
  | awk -F "_noControl" '{print $1}' \
  | parallel -j 8 "FRiP_ChIP {}_noControl_peaks.narrowPeak {}_dedup.bam ${MODE} 2>> {}.log" > FRiP_scores.txt
  
  fi
  
(>&2 paste -d " " <(echo '[INFO] FRiP score calculation ended on') <(date))

####################################################################################################################################

## Summary report:
(>&2 paste -d " " <(echo '[INFO] MultiQC started on') <(date))

multiqc -o multiqc_all ./ 

(>&2 paste -d " " <(echo '[INFO] MultiQC ended on') <(date))
