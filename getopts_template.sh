#!/bin/bash

## A template for getopts as used in the DNAseq_lowlevel script,
## thanks to jrjhealey for the template
## see https://github.com/jrjhealey/BadassBash/blob/master/BadassBash.sh#L50

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

-h | --help     : Show this message
-g | --genome   : The reference genome                   [mm10, hg38]
-a | --atac-seq : turn on filtering options for ATAC-seq [atacseq, chipseq]
-j | --jobs     : Number of parallel jobs                [4]
-i | --input    : The input file type                    [fq_single, fq_paired, 
                                                         bam_single, bam_paired,
                                                         cram_single, cram_paired]
-t | --threads  : Number of threads per job              [16]
-d | --alnOnly  : No peak calling / bigwig creation
-m | --macs     : Path to the macs2 executable
-r | --rscript  : Path to Rscript executable
-n | --genrich  : Whether to call peaks with Genrich
                 for ATAC-seq data
-c | --macsPeak : Whether to call peaks with macs
                 to get an idea of data quality

---------------------------------------------------------------------------------------------
'
}; if [[ -z "$1" ]] || [[ $1 == -h ]] || [[ $1 == --help ]]; then usage; exit; fi
	
#########################################################################################################

for arg in "$@"; do                         # for every arg in the commandline array ("$@")
 shift                                      # Shift by one so as to skip the script name
 case "$arg" in
   "--help")        set -- "$@" "-h"   ;;   # 
   "--genome")      set -- "$@" "-g"   ;;   # 
   "--assay")       set -- "$@" "-a"   ;;   #
   "--jobs")        set -- "$@" "-j"   ;;   #
   "--input")       set -- "$@" "-i"   ;;   #
   "--threads")     set -- "$@" "-t"   ;;   #
   "--alnOnly")     set -- "$@" "-d"   ;;   #
   "--macs")        set -- "$@" "-m"   ;;   #
   "--rscript")     set -- "$@" "-r"   ;;   #
   "--genrich")     set -- "$@" "-n"   ;;   #
   "--macsPeak")    set -- "$@" "-c"   ;;   #
   *)             set -- "$@" "$arg" ;;     # Lastly, deal with any unmatched args.
 esac
done
	
## Call getopts to assign the arguments to variables for use in the script:

while getopts hncmrd:g:a:j:i:t: OPT     # args followed by ":" expect an argument
  do   
  case ${OPT} in                           # Letters followed by a ":" are expecting an argument.
    g) Genome="${OPTARG}"    ;;
    a) Assay="${OPTARG}"     ;;
    j) Jobs="TRUE"           ;;
    d) AlnOnly="${OPTARG}"   ;;
    i) Inputtype="${OPTARG}" ;;
    t) Threads="${OPTARG}"   ;;
    m) MACS2="${OPTARG}"     ;;
    r) RSCRIPT="${OPTARG}"   ;;
    n) DoGenrich="TRUE" ;;
    c) DoMacs="TRUE"    ;;
  esac
done	

#########################################################################################################

## Check for args, set defaults if empty:
if [[ "${Jobs}" == "" ]]; then Jobs=16; fi

if [[ "${Threads}" == "" ]]; then Threads=14; fi

if [[ "${AlnOnly}" != "TRUE" ]]; then AlnOnly="FALSE"; fi

if [[ "${MACS2}" == "" ]]; then
  MACS2="$HOME/anaconda3_things/anaconda3/envs/py27_env/bin/macs2"
  fi
  
if [[ "${RSCRIPT}" == "" ]]; then
  RSCRIPT="$HOME/anaconda3_things/anaconda3/envs/R_env/bin/Rscript"
  fi
  
if [[ "${DoGenrich}" != "TRUE" ]]; then
  DoGenrich="FALSE"
  fi
  
if [[ "${DoMacs}" != "TRUE" ]]; then
  DoMacs="FALSE"
  fi  
  
#########################################################################################################

## Return summary:

echo ''
echo '---------------------------------------------------------------------------------------------'
echo '[Info] Running with these parameters:'
echo '       --genome   = '"${Genome}"
echo '       --assay    = '"${Assay}"
echo '       --jobs     = '"${Jobs}"
echo '       --input    = '"${Inputtype}"
echo '       --threads  = '"${Threads}"
echo '       --macs     = '"${MACS2}"
echo '       --rscript  = '"${RSCRIPT}"
echo '       --genrich  = '"${DoGenrich}"
echo '       --macsPeak = '"${DoMacs}"
echo '       --alnOnly  = '"${AlnOnly}"
echo '---------------------------------------------------------------------------------------------'
echo ''

#########################################################################################################
