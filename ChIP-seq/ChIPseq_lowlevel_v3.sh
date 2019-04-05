#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --partition=normal
#SBATCH --mem=80G
#SBATCH --time=24:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=ChIPseq_lowlevel

######################################################################################################################################

## Script for lowlevel processing of ChIP-seq data: Alignment, filtering, quality control, browser tracks:

######################################################################################################################################

## Due to different conda envs (macs2 needs py27, and Rscript is in R env), some things have to be set explicitely:

GENOME="mm10"
MODE="SE"
BLACKLIST="/scratch/tmp/a_toen03/Genomes/mm10/mm10_consensusBL.bed"
RSCRIPT="$HOME/anaconda3_things/anaconda3/envs/R_env/bin/Rscript"
MACS="$HOME/anaconda3_things/anaconda3/envs/py27_env/bin/macs2"

if [[ ${GENOME} == "hg38" ]]; then
	IDX="/scratch/tmp/a_toen03/Genomes/hg38/bowtie2_index_noALT_withDecoy/hg38_noALT_withDecoy.fa"
	fi
  
if [[ ${GENOME} == "mm10" ]]; then
	IDX="/scratch/tmp/a_toen03/Genomes/mm10/bowtie2_idx/mm10"
	fi  

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
       bigWigToBedGraph $MACS bg2bw pigz bamcollate2 idr)

for i in $(echo ${TOOLS[*]}); do
  PathCheck $i; done
  
if [[ -e missing_tools.txt ]] && [[ $(cat missing_tools.txt | wc -l | xargs) > 0 ]]; then
  echo '[ERROR] Tools missing in PATH -- see missing_tools.txt' && exit 1
  fi

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

function Fq2Bam {

  BASENAME=$1
  TYPE=$2
  IDX=$3
	
  if [[ $TYPE == "PE" ]]; then
    IND="paired-end mode"
    fi
  if [[ $TYPE == "SE" ]]; then
    IND="single-end mode"
    fi
			
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'in' ${IND} 'started on') <(date))
  
  if [[ $TYPE == "PE" ]]; then
    if [[ ! -e ${BASENAME}_1.fastq.gz ]] || [[ ! -e ${BASENAME}_2.fastq.gz ]]; then
    echo '[ERROR] At least one input files is missing -- exiting' && exit 1
      fi
    else
      if [[ ! -e ${BASENAME}.fastq.gz ]] ; then
        echo '[ERROR] At least one input files is missing -- exiting' && exit 1
        fi
    fi   
    
  ## Nextera adapter:
  ADAPTER="AGATCGGAAGAGC"
    
  PAIREDRUN="seqtk mergepe ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz | \
  	     cutadapt -j 1 -a ${ADAPTER} -A ${ADAPTER} --interleaved -m 18 --max-n 0.1 - | \
	     bowtie2 --very-sensitive --threads 16 -x $IDX --interleaved - | \
	     samtools fixmate -m -@ 2 -O SAM - -"
	
  SINGLERUN="cutadapt -j 1 -a ${ADAPTER} -m 18 --max-n 0.1 ${BASENAME}.fastq.gz | \
	     bowtie2 --very-sensitive --threads 16 -x $IDX -U -"
	
  if [[ $TYPE == "PE" ]]; then
    RUN=${PAIREDRUN}
      fi
  if [[ $TYPE == "SE" ]]; then
    RUN=${SINGLERUN}
    fi
	
  ## ALign and filter:
  eval ${RUN} | \
  samblaster --ignoreUnmated | \
  sambamba view -f bam -S -l 5 -t 2 -o /dev/stdout /dev/stdin | \
  tee >(sambamba flagstat -t 2 /dev/stdin > ${BASENAME}_raw.flagstat) | \
  sambamba sort -m 5G --tmpdir=./ -l 5 -t 16 -o ${BASENAME}_raw.bam /dev/stdin  
  
  BamCheck ${BASENAME}_raw.bam
    
  ## 1) remove non-primary chromosomes, low qual. and non-primary alignments, but keep duplicates:
  if [[ ! -e tmp_chromSizes.txt ]]; then
    samtools idxstats ${BASENAME}_raw.bam > tmp_chromSizes.txt
    fi
  
  if [[ $TYPE == "PE" ]]; then
    NUMFILTER=1
    fi
  if [[ $TYPE == "SE" ]]; then
    NUMFILTER=0
    fi
		
  samtools idxstats ${BASENAME}_raw.bam | cut -f 1 | grep -vE 'chrM|_random|chrU|chrEBV|\*' | \
  xargs sambamba view -l 5 -f bam -t 8 --num-filter=${NUMFILTER}/2308 --filter='mapping_quality > 19' \
  -o /dev/stdout ${BASENAME}_raw.bam | \
  tee ${BASENAME}_dup.bam | \
  tee >(samtools index - ${BASENAME}_dup.bam.bai) | \
  sambamba view -l 5 -f bam -t 8 --num-filter=/1024 -o ${BASENAME}_dedup.bam /dev/stdin | \
  
  BamCheck ${BASENAME}_dup.bam
  
  ls ${BASENAME}*dup.bam | parallel "sambamba flagstat -t 8 {} > {.}.flagstat"
  
  BamCheck ${BASENAME}_dedup.bam
  
  bamCoverage --bam ${BASENAME}_dedup.bam -p 16 -o ${BASENAME}_dedup_CPM.bigwig --normalizeUsing CPM -bs 1 -e 200
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' $1 'in' ${IND} 'ended on') <(date))

}

export -f Fq2Bam

######################################################################################################################################

## Count fraction of cutting sites in Genrich peaks:
function FRiP {
  
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

}; export -f FRiP

####################################################################################################################################

## fastqc:
ls *fastq.gz | parallel -j 70 "fastqc -t 1 {}"

####################################################################################################################################

## Aln/filtering:
if [[ $MODE != "PE" ]] && [[ $MODE != "SE" ]]; then
  echo '[ERROR] Missing SE/PE information in Fq2Bam function -- exiting'
  exit 1; fi
 
if [[ $MODE == "PE" ]]; then
  ls *_1.fastq.gz | awk -F "_1.fastq.gz" '{print $1}' | sort -u | parallel -j 4 "Fq2Bam {} ${MODE} ${IDX} 2>> {}.log"
  fi

if [[ $MODE == "SE" ]]; then
  ls *.fastq.gz | awk -F ".fastq.gz" '{print $1}' | sort -u | parallel -j 4 "Fq2Bam {} ${MODE} ${IDX} 2>> {}.log"
  fi

####################################################################################################################################

## Insert Sizes given paired-end data:
if [[ $MODE == "PE" ]]; then
  (>&2 paste -d " " <(echo '[INFO] CollectInsertSizes started on') <(date))

  ls *_dedup.bam | \
    parallel "picard CollectInsertSizeMetrics I={} O={.}_InsertSizes.txt H={.}_InsertSizes.pdf QUIET=true VERBOSITY=ERROR 2> /dev/null"
  
  (>&2 paste -d " " <(echo '[INFO] CollectInsertSizes started on') <(date))
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

## Peaks for each sample at 1% FDR:
if [[ $GENOME == "mm10" ]]; then GFLAG="mm"; fi
if [[ $GENOME == "hg38" ]]; then GFLAG="hs"; fi

## Call Peaks:
if [[ $MODE == "PE" ]]; then FORMAT=BAMPE; fi
if [[ $MODE == "SE" ]]; then FORMAT=BAM; fi

ls *_dedup.bam | awk -F "_dedup.bam" '{print $1}' | sort -u | \
  parallel "$MACS callpeak -t {}_dedup.bam -n {}_noControl -g $GFLAG -f ${FORMAT}" 2> macs_reports.log  
		  
####################################################################################################################################

## FRiP:
(>&2 paste -d " " <(echo '[INFO] FRiP score calculation started on') <(date))

ls *_noControl_peaks.narrowPeak | \
  awk -F "_noControl" '{print $1}' | \
  parallel -j 4 "FRiP {}_noControl_peaks.narrowPeak {}_dedup.bam ${MODE} 2>> {}.log" > FRiP_scores.txt
  
(>&2 paste -d " " <(echo '[INFO] FRiP score calculation endedn on') <(date))

####################################################################################################################################

## Summary report:
(>&2 paste -d " " <(echo '[INFO] MultiQC started on') <(date))
multiqc -o multiqc_all ./ 
(>&2 paste -d " " <(echo '[INFO] MultiQC ended on') <(date))
