#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --partition=normal
#SBATCH --mem=80G
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=HiC_alignments

######################################################################################################################################
######################################################################################################################################

## Alignment of Hi-C data

######################################################################################################################################
######################################################################################################################################

## => JOBS can be set to 4 if using the SMP nodes or 2 if on the normal Skylake nodes

GENOME="mm10"
JOBS=2

if [[ ${GENOME} == "hg38" ]]; then
  IDX="/scratch/tmp/a_toen03/Genomes/hg38/bowtie2_index_noALT_withDecoy/hg38_noALT_withDecoy.fa"
  fi
  
if [[ ${GENOME} == "mm10" ]]; then
  IDX="/scratch/tmp/a_toen03/Genomes/mm10/bowtie2_idx/mm10"
  fi  
  
######################################################################################################################################
######################################################################################################################################

## => Section to check parameters and tools being in PATH:

if [[ ! $(echo $GENOME | grep -E 'hg38|mm10') ]]; then 
  echo '[ERROR] GENOME is neither of the supported hg38 and mm10 -- exiting' && exit 1
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
TOOLS=(bowtie2 cutadapt samtools multiqc fastqc)

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

function Fq2Bam {
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' "${1}" 'started on') <(date))
  
  BASENAME="${1}"
  IDX="${2}"

  #############################################################################################################
  
  ## Align each mate pair with bowtie2:
  ADAPTER="AGATCGGAAGAGC"
  
  cutadapt -a "${ADAPTER}" -m 18 --max-n 0.1 "${BASENAME}".fastq.gz \
    | bowtie2 --reorder --threads 32 -x "${IDX}" -U - \
    | samblaster --ignoreUnmated \
    | tee >(samtools flagstat - > "${BASENAME}".flagstat) \
    | samtools view -@ 4 -o "${BASENAME}".bam

  ## Check if alignment went ok based on BAM file integrity:
  BamCheck "${BASENAME}".bam
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' "${1}" 'ended on') <(date))

}

export -f Fq2Bam

######################################################################################################################################

## Run:
ls *.fastq.gz \
  | awk -F ".fastq.gz" '{print $1}' \
  | sort -u \
  | parallel -j $JOBS "Fq2Bam {} ${IDX} 2>> {}.log"
  
## fastqc:  
ls *.fastq.gz | parallel "fastqc {}"  

## multiqc summary:
multiqc .
