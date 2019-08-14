#!/bin/bash

## Help section:
usage(){
echo '

## DNA-seq lowlevel processing for ATACseq/ChIPseq-like data
 
   Naming conventions are:
   1) single-end fastq : Basename.fastq.gz
   2) paired-end fastq : Basename_1/2.fastq.gz
   3)            ubam  : Basename_ubam.bam
   4)            ucram : Basename_ucram.cram
   
## Options:

---------------------------------------------------------------------------------------------

-h | --help      : Show this message
-g | --genome    : The reference genome                   [mm10, hg38]
-a | --atacseq  : turn on filtering options for ATAC-seq 
-j | --jobs      : Number of parallel jobs                [4]
-s | --seqtype   : The input file type                    [fq_single, fq_paired, 
                                                           bam_single, bam_paired,
                                                           cram_single, cram_paired]
-t | --threads   : Number of threads per job              [16]
-d | --alnOnly   : Exit after alignment
-m | --macs      : Path to the MACS executable            [default the one on py27_env on HPC]
-c | --macsPeaks : Whether to call peaks with macs
                   to get an idea of data quality

---------------------------------------------------------------------------------------------
'
}; if [[ -z "$1" ]] || [[ $1 == -h ]] || [[ $1 == --help ]]; then usage; exit; fi
	
#################################################################################################################################################

## Thanks to jrjhealey for the template
## https://github.com/jrjhealey/BadassBash/blob/master/BadassBash.sh#L50

for arg in "$@"; do                         # for every arg in the commandline array ("$@")
 shift                                      # Shift by one so as to skip the script name
 case "$arg" in
   "--genome")      set -- "$@" "-g"   ;;   # 
   "--atacseq")     set -- "$@" "-a"   ;;   #
   "--jobs")        set -- "$@" "-j"   ;;   #
   "--seqtype")     set -- "$@" "-s"   ;;   #
   "--threads")     set -- "$@" "-t"   ;;   #
   "--macs")        set -- "$@" "-m"   ;;   #
   "--macsPeaks")   set -- "$@" "-c"   ;;   #
   "--alnOnly")     set -- "$@" "-d"   ;;   #
   *)               set -- "$@" "$arg" ;;   # Lastly, deal with any unmatched args.
 esac
done
	
## Call getopts to assign the arguments to variables for use in the script:

while getopts cmdag:j:s:t: OPT         # args followed by ":" expect an argument
  do   
  case ${OPT} in                            # Letters followed by a ":" are expecting an argument.
    g) Genome="${OPTARG}"    ;;
    a) ATACseq="TRUE"        ;;             # boolean
    j) Jobs="${OPTARG}"      ;;
    s) Inputtype="${OPTARG}" ;;
    t) Threads="${OPTARG}"   ;;
    m) MACS="${OPTARG}"      ;;
    c) DoMacs="TRUE"         ;;
    d) AlnOnly="TRUE"        ;;
  esac
done	

#################################################################################################################################################

## Check for args, set defaults if empty:
if [[ "${ATACseq}" != "TRUE" ]]; then ATACseq="FALSE"; fi

if [[ "${Jobs}" == "" ]]; then Jobs=4; fi

if [[ "${Threads}" == "" ]]; then Threads=16; fi

if [[ "${MACS}" == "" ]]; then
  MACS="$HOME/anaconda3_things/anaconda3/envs/py27_env/bin/macs2"; fi
  
if [[ ${DoMacs} != "TRUE" ]]; then DoMacs="FALSE"; fi

if [[ ${AlnOnly} != "TRUE" ]]; then AlnOnly="FALSE"; fi

if [[ ${AlnOnly} == "TRUE" ]]; then 
  DoMacs="FALSE"; fi
  
#################################################################################################################################################

## Return summary:

echo ''
echo '---------------------------------------------------------------------------------------------'
echo '[Info] Running with these parameters:'
echo '       --genome   = '"${Genome}"
echo '       --atacseq  = '"${ATACseq}"
echo '       --jobs     = '"${Jobs}"
echo '       --seqtype  = '"${Inputtype}"
echo '       --threads  = '"${Threads}"
echo '       --macs     = '"${MACS}"
echo '       --macsPeak = '"${DoMacs}"
echo '       --alnOnly  = '"${AlnOnly}"
echo '---------------------------------------------------------------------------------------------'
echo ''

#################################################################################################################################################
	
if [[ ${Genome} == "hg38" ]]; then
  Idx="/scratch/tmp/a_toen03/Genomes/hg38/bowtie2_index_noALT_withDecoy/hg38_noALT_withDecoy.fa"
  BLACKLIST="/scratch/tmp/a_toen03/Genomes/hg38/Blacklists/hg38_consensusBL.bed"
  fi
  
if [[ ${Genome} == "mm10" ]]; then
  Idx="/scratch/tmp/a_toen03/Genomes/mm10/bowtie2_idx/mm10"
  BLACKLIST="/scratch/tmp/a_toen03/Genomes/mm10/mm10_consensusBL.bed"
  fi  
  
#################################################################################################################################################

## Check if required tools are in PATH and/or callable:

if [[ -e missing_tools.txt ]]; then rm missing_tools.txt; fi

## Function that checks if required tools are callable:
function PathCheck {
  
  if [[ $(command -v $1 | wc -l | xargs) == 0 ]]; then 
    echo ${1} >> missing_tools.txt
  fi
    
}; export -f PathCheck

## All tools:
TOOLS=(seqtk cutadapt bwa samtools samblaster sambamba bedtools mawk bgzip tabix \
       featureCounts bc bamCoverage parallel fastqc picard multiqc \
       bigWigToBedGraph ${MACS} bg2bw pigz featureCounts)

## Loop through the tools and write to missing_tools.txt all those not in PATH / callable:
for i in $(echo ${TOOLS[*]}); do
  PathCheck $i; done
  
## If tools are missing write them to <missing_tools.txt>
if [[ -e missing_tools.txt ]] && [[ $(cat missing_tools.txt | wc -l | xargs) > 0 ]]; then
  echo '[ERROR] Tools missing in PATH -- see missing_tools.txt' && exit 1; fi

#################################################################################################################################################

## Exit function if BAM file looks corrupted or is missing after a step:
function ExitBam {
  (>&2 echo '[ERROR]' "${1}" 'looks suspicious or is empty -- exiting') && exit 1
}; export -f ExitBam

#################################################################################################################################################

## Check if BAM is corrupted or empty:
function BamCheck {
  
  Basename="${1%_*}"
  samtools quickcheck -q $1 && echo '' >/dev/null || ExitBam "${1}"
    
  ## Also check if file is not empty:
  if [[ $(samtools view "${1}" | head -n 1 | wc -l) < 1 ]]; then
    ExitBam $Basename; fi
  
}; export -f BamCheck  

#################################################################################################################################################

## Function to get the % of reads mapped to chrM:
function mtDNA {

  BamCheck $1
    
  if [[ ! -e "${1}".bai ]];
    then sambamba index -t 8 ${1}; fi
      
  mtReads=$(samtools idxstats $1 | grep 'chrM' | cut -f 3)
  totalReads=$(samtools idxstats $1 | awk '{SUM += $3} END {print SUM}')

  echo '[mtDNA Content]' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%' > "${1%.bam}"_mtDNA.txt
  
}; export -f mtDNA

#################################################################################################################################################

## Alignment /filtering function:

function Fq2Bam {

  Basename="${1}"
  Inputtype="${2}"
  Idx="${3}"
  ATACseq="${4}"
  Threads="${5}"
	
#################################################################################################################################################
  
  ## default adapter types
  if [[ ${ATACseq} == "TRUE" ]]
    then
      ADAPTER="CTGTCTCTTATACACATCT"
    else
      ADAPTER="AGATCGGAAGAGC"
    fi  
  
  ## six different commands depending on the sequencing type (se/pe) and the input (fq/bam/cram)
  CUT_FQSE="cutadapt -j 1 -a "${ADAPTER}" -m 18 --max-n 0.1 "${Basename}".fastq.gz"
	
  CUT_FQPE="seqtk mergepe "${Basename}"_1.fastq.gz "${Basename}"_2.fastq.gz \
  	        | cutadapt -j 1 -a "${ADAPTER}" -A "${ADAPTER}" --interleaved -m 18 --max-n 0.1 -"
  	        
  CUT_SBAM="samtools fastq -n -@ 2 "${Basename}"_ubam.bam \
            | cutadapt -j 1 -a "${ADAPTER}" -m 18 --max-n 0.1 -"	 
            
  CUT_PBAM="samtools fastq -n -@ 2 "${Basename}"_ubam.bam \
            | cutadapt -j 1 -a "${ADAPTER}" -A "${ADAPTER}" --interleaved -m 18 --max-n 0.1 -"
            
  CUT_SCR="samtools fastq -n -@ 2 "${Basename}"_ucram.bam \
            | cutadapt -j 1 -a "${ADAPTER}" -m 18 --max-n 0.1 -"	 
            
  CUT_PCR="samtools fastq -n -@ 2 "${Basename}"_ucram.cram \
            | cutadapt -j 1 -a "${ADAPTER}" -A "${ADAPTER}" --interleaved -m 18 --max-n 0.1 -"          
  
#################################################################################################################################################
  
  ## Bowtie2 basic commands depending on $TYPE, reading from stdin:
  
  BOWTIE2PE="bowtie2 --very-sensitive --threads ${Threads} -X 2000 --rg-id "${Basename}" -x "${Idx}" --interleaved - \
             | samtools fixmate -m -@ 2 -O SAM - -"

  BOWTIE2SE="bowtie2 --very-sensitive --threads ${Threads} --rg-id "${Basename}" -x "${Idx}" -U -"             
	
#################################################################################################################################################
  
  ## Check mode and set concat command accordingly:
  
  if [[ ${Inputtype} == "fq_single" ]]; then COMMAND1="${CUT_FQSE} | ${BOWTIE2SE}"; fi
  
  if [[ ${Inputtype} == "fq_paired" ]]; then COMMAND1="${CUT_FQPE} | ${BOWTIE2PE}"; fi 
  
  if [[ ${Inputtype} == "bam_single" ]]; then COMMAND1="${CUT_SBAM} | ${BOWTIE2SE}"; fi
  
  if [[ ${Inputtype} == "bam_paired" ]]; then COMMAND1="${CUT_PBAM} | ${BOWTIE2PE}"; fi
  
  if [[ ${Inputtype} == "cram_single" ]]; then COMMAND1="${CUT_SBAM} | ${BOWTIE2SE}"; fi
  
  if [[ ${Inputtype} == "cram_paired" ]]; then COMMAND1="${CUT_PBAM} | ${BOWTIE2PE}"; fi
  
  ## Summary message:    
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' "${1}" 'in' ${Inputtype} 'mode started on') <(date))

  ## Run the lowlevel alignment pipeline:
  eval "${COMMAND1}" \
  | samblaster --ignoreUnmated \
  | sambamba view -f bam -S -l 5 -t 2 -o /dev/stdout /dev/stdin \
  | tee >(sambamba flagstat -t 2 /dev/stdin > "${Basename}"_raw.flagstat) \
  | sambamba sort -m 5G --tmpdir=./ -l 5 -t 16 -o "${Basename}"_raw.bam /dev/stdin  
  
  ## Check if alignment went ok based on BAM file integrity:
  BamCheck "${Basename}"_raw.bam
  
  ## Get mtDNA percentage if ATAC-seq:
  if [[ ${ATACseq} == "TRUE" ]]; then mtDNA "${Basename}"_raw.bam; fi
  
  ## Start filtering the BAM file:
  if [[ ! -e tmp_chromSizes.txt ]]; then
    samtools idxstats "${Basename}"_raw.bam > tmp_chromSizes.txt
    fi
  
  ## Decide appropriate flag to keep only mapped reads depending on paired/single-end sequencing:
  if [[ $(echo ${Inputtype} | grep '_paired') ]]; then
    NUMFILTER=1; fi
    
  if [[ $(echo ${Inputtype} | grep '_single') ]]; then
    NUMFILTER=0; fi
    
	## Now define filtering that happens for all DNA-seq:
  GENERALFILTER="samtools idxstats "${Basename}"_raw.bam \
                 | cut -f 1 \
                 | grep -vE 'chrM|_random|chrU|chrEBV|\*' \
                 | xargs sambamba view -l 5 -f bam -t 8 --num-filter="${NUMFILTER}"/2308 --filter='mapping_quality > 19' \
                   -o /dev/stdout "${Basename}"_raw.bam \
                 | tee "${Basename}"_dup.bam \
                 | tee >(samtools index - "${Basename}"_dup.bam.bai)"
     
  ## Cat commands together:            	
  if [[ ${ATACseq} == "FALSE" ]]; then
    eval "${GENERALFILTER}" \
    | sambamba view -l 5 -f bam -t 8 --num-filter=/1024 -o ${Basename}_dedup.bam /dev/stdin
    fi
  
  if [[ ${ATACseq} == "TRUE" ]]; then
  
    eval "${GENERALFILTER}" \
    | sambamba view -l 5 -f bam -t 8 --num-filter=/1024 -o /dev/stdout /dev/stdin \
    | tee >(tee "${Basename}"_dedup.bam | samtools index - "${Basename}"_dedup.bam.bai) \
    | bedtools bamtobed -i - \
    | mawk 'OFS="\t" {if ($6 == "+") print $1, $2+4, $2+5, ".", ".", $6} {if ($6 == "-") print $1, $3-5, $3-4, ".", ".", $6}' \
    | sort -k1,1 -k2,2n -k3,3n -k6,6 -S10G --parallel=10 \
    | tee >(bgzip -@ 6 > "${Basename}"_cutsites.bed.gz) \
    | bedtools genomecov -bg -i - -g tmp_chromSizes.txt \
    | bg2bw -i /dev/stdin -c tmp_chromSizes.txt -o "${Basename}"_cutsites.bigwig
    
  fi
  
  BamCheck "${Basename}"_dup.bam
  BamCheck "${Basename}"_dedup.bam
    
  if [[ $(echo ${Inputtype} | grep '_paired') ]]; then
    picard CollectInsertSizeMetrics \
    I="${Basename}"_dedup.bam \
    O="${Basename}"_InsertSizes.txt \
    H="${Basename}"_InsertSizes.pdf \
    QUIET=true VERBOSITY=ERROR VALIDATION_STRINGENCY=LENIENT 2> /dev/null
  fi
  
  ## Flagstat that I couldn't squeeze into the tee part above:
  tmpthread=$(bc <<< "${Threads}/2")
  if [[ ${tmpthread} < 1 ]]; then tmpthread=1; fi
  ls "${Basename}"*dup.bam | parallel -j 2 "sambamba flagstat -t ${tmpthread} {} > {.}.flagstat"
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' "${1}" 'in' ${Inputtype} 'mode ended on') <(date))

}

export -f Fq2Bam

#################################################################################################################################################

##

#################################################################################################################################################

## Run Fq2Bam:

if [[ "${Inputtype}" == "fq_single" ]]; then 
  FileDefine=$(ls *.fastq.gz 2> /dev/null | awk -F ".fastq.gz" '{print $1 | "sort -u"}')
fi

if [[ "${Inputtype}" == "fq_paired" ]]; then 
  FileDefine=$(ls *_1.fastq.gz 2> /dev/null | awk -F "_1.fastq.gz" '{print $1 | "sort -u"}')
fi
  
if [[ "${Inputtype}" == "bam_single" ]] || [[ "${Inputtype}" == "bam_paired" ]] ; then 
  FileDefine=$(ls *_ubam.bam 2> /dev/null | awk -F "_ubam.bam" '{print $1 | "sort -u"}')
fi  
  
if [[ "${Inputtype}" == "cram_single" ]] || [[ "${Inputtype}" == "cram_paired" ]] ; then 
  FileDefine=$(ls *_ucram.cram 2> /dev/null | awk -F "_ucram.cram" '{print $1 | "sort -u"}')
fi
  
if [[ $(echo $FileDefine | wc -l | xargs) == 0 ]]; then 
  echo "[Error]: No input files found!"
  exit 1
fi  
  
## Start the job:  
echo "${FileDefine}" \
  | tr " " "\n" \
  | parallel -j "${Jobs}" "Fq2Bam {} ${Inputtype} ${Idx} ${ATACseq} ${Threads} 2> {}.log"

#################################################################################################################################################

## fastqc of raw output:
echo "${FileDefine}" \
  | tr " " "\n" \
  | parallel -j "$(nproc)" "fastqc {}_raw.bam"

#################################################################################################################################################

## Call peaks with macs:

if [[ ${DoMacs} == "TRUE" ]]; then
  
  (>&2 paste -d " " <(echo '[INFO] Peaks/FRiPs started on') <(date))
  
  ## effective genome size for macs:
  if [[ "${Genome}" == "mm10" ]]; then GFLAG="mm"; fi
  if [[ "${Genome}" == "hg38" ]]; then GFLAG="hs"; fi
  
  ## MACS for ATAC-seq using cutting sites at 1% FDR:
  if [[ "${ATACseq}" == "TRUE" ]]; then
  
    if [[ $(ls *_cutsites.bed.gz 2> /dev/null | wc -l | xargs) == 0 ]]; then
      echo "[Error]: No input files found for peak calling"
      exit 1
    fi
    
    ls *_cutsites.bed.gz \
      | awk -F "_cutsites.bed.gz" '{print $1 | "sort -u"}' \
      | parallel "$MACS callpeak -t {}_cutsites.bed.gz -n {} -g $GFLAG \
                  --extsize 100 --shift -50 --nomodel --keep-dup=all -f BED --call-summits -q 0.01 \
    		          --min-length 150 2>> {}.log"
    		          
  fi
  
  ## ChIP-seq with default settings:
  if [[ "${ATACseq}" == "FALSE" ]]; then
  
    if [[ $(ls *_dedup.bam 2> /dev/null | wc -l | xargs) == 0 ]]; then
      echo "[Error]: No input files found for peak calling"
      exit 1
    fi
    
    if [[ $(echo $Inputtype | grep '_paired') ]]; then
      BAM="BAMPE"
    else
      BAM="BAM"
    fi  
  
    ls *_dedup.bam \
      | awk -F "_dedup.bam" '{print $1 | "sort -u"}' \
      | parallel "$MACS callpeak -t {}_dedup.bam -n {} -g $GFLAG --keep-dup=all -f ${BAM} 2>> {}.log"
              
  fi
  
  ## move things to special directory
  if [[ ! -e ./MacsDir ]]; then mkdir ./MacsDir; fi
    
    ls \
      | grep "_peaks.narrowPeak" \
      | while read p; do
          bedtools intersect -v -a ${p} -b ${BLACKLIST} \
          | tee ./MacsDir/${p} \
          | awk 'OFS="\t" {print $1"_"$2+1"_"$3, $1, $2+1, $3, "+"}' > ${p%.narrowPeak}.saf
        done < /dev/stdin
    rm *_summits.bed *.xls *.narrowPeak
  
  ## FRiPs:
  if [[ $(echo ${Inputtype} | grep '_paired') ]]; then
    P="-p"
  else
      P=""
  fi
  
  ls *.saf \
    | awk -F "_peaks.saf" '{print $1 | "sort -u"}' \
    | parallel -j ${Jobs} "featureCounts ${P} -a {}_peaks.saf -F SAF -T ${Threads} -o {}_countmatrix.txt {}_dedup.bam 2>> {}.log"
  
  ## calcualte frips from featureCounts output:
  for b in *_countmatrix.txt; do
    
    ASSIGNED=$(grep -w 'Assigned' ${b}.summary | cut -f2)
    UNASSIGNED=$(grep -w 'Unassigned_NoFeatures' ${b}.summary | cut -f2)
    
    paste <(echo "${b%_countmatrix.txt}") <(bc <<< "scale=6;${ASSIGNED}/(${ASSIGNED}+${UNASSIGNED})") >> FRiPs.txt
    
  done

  (>&2 paste -d " " <(echo '[INFO] Peaks/FRiPs ended on') <(date))
  
fi
  
#################################################################################################################################################  
  
## Summary report:
multiqc -o multiqc_all ./ 
