#!/bin/bash

## Wrapper to quantify fastq files (RNA-seq) with salmon (v.1.0.0 and later / pufferfish index)

## Help section:
usage(){
echo '

## RNA-seq quantification with salmon v.1.0.0 and later:
 
   Naming conventions are:
   1) single-end fastq : Basename.fastq.gz
   2) paired-end fastq : Basename_1/2.fastq.gz

## Options:

---------------------------------------------------------------------------------------------

-h | --help      : Show this message
-i | --idx       : the transcriptome index
-m | --mode      : single or paired-end data                       [single,paired]
-n | --noLength  : turn off length correction if 3`-enriched
                   libraries are quantified                       
-t | --threads   : number of threads per sample                    [8]
-j | --jobs      : number of parallel jobs                         [8]
-l | --libtype   : library type                                    [A]
-f | --fastqc    : guess what this does                            [TRUE]
-u | --trim      : use cutadapt to trim adapters                   [FALSE]
---------------------------------------------------------------------------------------------
'
}; if [[ -z "$1" ]] || [[ $1 == -h ]] || [[ $1 == --help ]]; then usage; exit; fi
	
#################################################################################################################################################

## Thanks to jrjhealey for the template
## https://github.com/jrjhealey/BadassBash/blob/master/BadassBash.sh#L50

for arg in "$@"; do                         # for every arg in the commandline array ("$@")
 shift                                      # Shift by one so as to skip the script name
 case "$arg" in
   "--idx")         set -- "$@" "-i"   ;;   # 
   "--mode")        set -- "$@" "-m"   ;;   # 
   "--noLength")    set -- "$@" "-n"   ;;   #
   "--jobs")        set -- "$@" "-j"   ;;   #
   "--threads")     set -- "$@" "-t"   ;;   #
   "--libtype")     set -- "$@" "-l"   ;;   #
   "--fastqc")      set -- "$@" "-f"   ;;   #
   "--trim")        set -- "$@" "-u"   ;;   #
   *)               set -- "$@" "$arg" ;;   # Lastly, deal with any unmatched args.
 esac
done
	
## Set default and call getopts to assign the arguments to variables for use in the script:
NoLength="FALSE"
Trim="FALSE"
Fastqc="FALSE"

while getopts i:m:j:t:l:nfu OPT      # args followed by ":" expect an argument
  do   
  case ${OPT} in                     # Letters followed by a ":" are expecting an argument.
    i) Idx="${OPTARG}"       ;;
    m) Mode="${OPTARG}"      ;;
    n) NoLength="TRUE"       ;;     
    j) Jobs="${OPTARG}"      ;;
    t) Threads="${OPTARG}"   ;;
    l) LibType="${OPTARG}"   ;;
    f) Fastqc="TRUE"         ;;
    u) Trim="TRUE"           ;;
  esac
done	

#################################################################################################################################################

## Check for args, set defaults if empty:

if [[ "${Threads}" == "" ]]; then Threads=8; fi

if [[ "${Jobs}" == "" ]]; then Jobs=8; fi

if [[ "${LibType}" == "" ]]; then LibType=A; fi

if [[ "${Trim}" == "TRUE" ]]; then
  echo 'Trimming not yet implemented, do it manually before running this wrapper!'
  exit 1
fi

OPTIONS=(Idx Mode NoLength Jobs Threads LibType Fastqc Trim)
for i in $(echo ${OPTIONS[*]}); do
  export $i
done  

#################################################################################################################################################

## Return summary:

echo ''
echo '----------------------------------------------------------------------------------------------------------------------------------'
echo '[Info] Running with these parameters:'
echo '       --idx      = '"${Idx}"
echo '       --mode     = '"${Mode}"
echo '       --noLength = '"${NoLength}"
echo '       --jobs     = '"${Jobs}"
echo '       --threads  = '"${Threads}"
echo '       --libtype  = '"${LibType}"
echo '       --fastqc   = '"${Fastqc}"
echo '       --trim     = '"${Trim}"
echo '----------------------------------------------------------------------------------------------------------------------------------'
echo ''

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
TOOLS=(salmon fastqc multiqc)

## Loop through the tools and write to missing_tools.txt all those not in PATH / callable:
for i in $(echo ${TOOLS[*]}); do
  PathCheck $i; done
  
## If tools are missing write them to <missing_tools.txt>
if [[ -e missing_tools.txt ]] && [[ $(cat missing_tools.txt | wc -l | xargs) > 0 ]]; then
  echo '[ERROR] Tools missing in PATH -- see missing_tools.txt' && exit 1; fi

#################################################################################################################################################

function SALMON {
  
  SalmonBasic="salmon quant -l ${LibType} -i ${Idx} -p ${Threads} --validateMappings --writeUnmappedNames --no-version-check"
  
  (>&2 paste -d " " <(echo '[INFO]' 'Running salmon for' "${1}" 'in' ${Mode} 'mode -- started on') <(date))
  
  ## Paired-End run:              
  if [[ ${Mode} == "paired" ]]; then
  
    if [[ ! -e ${1}_1.fastq.gz || ! -e ${1}_2.fastq.gz ]]; then
    echo '[ERROR]: At least on of the input files is missing for' $1 && exit 1
    fi
    
    ## standard full-length RNA-seq:
    if [[ "${NoLength}" == "FALSE" ]]; then
      eval "${SalmonBasic}" --seqBias --gcBias -o ${1}_salmon -1 ${1}_1.fastq.gz -2 ${1}_2.fastq.gz
    fi  
    
    ## 3'-enriched libraries so no length correction:
    if [[ "${NoLength}" == "TRUE" ]]; then
      eval "${SalmonBasic}" --gcBias --noLengthCorrection -o ${1}_salmon -1 ${1}_1.fastq.gz -2 ${1}_2.fastq.gz
    fi
  
  fi
  
  ## Single-end              
  if [[ ${Mode} == "single" ]]; then
    
    if [[ ! -e ${1}.fastq.gz ]]; then
    echo '[ERROR]: Input file is missing for' $1 && exit 1
    fi
    
    ## standard full-length RNA-seq:
    if [[ "${NoLength}" == "FALSE" ]]; then
      eval "${SalmonBasic}" --seqBias -o ${1}_salmon -r ${1}.fastq.gz
    fi  
    
    ## 3'-enriched libraries so no length correction:
    if [[ "${NoLength}" == "TRUE" ]]; then
      eval "${SalmonBasic}" --noLengthCorrection -o ${1}_salmon -r ${1}.fastq.gz
    fi 
  
  fi
  
  (>&2 paste -d " " <(echo '[INFO]' 'Running salmon for' "${1}" 'in' ${Mode} 'mode -- ended on') <(date))
 
}; export -f SALMON

if [[ ${Mode} == "paired" ]]; then 
  ls *_1.fastq.gz \
  | awk -F "_1.fastq.gz" '{print $1}' \
  | parallel -j ${Jobs} "SALMON {} 2> {}_salmon.log"
fi

if [[ ${Mode} == "single" ]]; then 
  ls *.fastq.gz \
  | awk -F ".fastq.gz" '{print $1}' \
  | parallel -j ${Jobs} "SALMON {} 2> {}_salmon.log"
fi

if [[ ${Fastqc} == "TRUE" ]]; then ls *.fastq.gz | parallel -j 64 "fastqc {}"; fi

multiqc .
