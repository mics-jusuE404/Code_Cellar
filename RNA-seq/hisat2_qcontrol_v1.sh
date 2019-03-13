#!/bin/bash

## Quality control and alignment of PE RNA-seq data with hisat2, RSeQC:

#######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --partition=normal
#SBATCH --mem=80G
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=hisat2_qualcontrol
#######

########################################################################################################################################

GENOME="mm10"
MODE="PE"

if [[ ${GENOME} == "mm10" ]]; then
  HISAT2_IDX="/scratch/tmp/a_toen03/Genomes/mm10/hisat2_index/mm10"
  SPLICE_FILE="/scratch/tmp/a_toen03/Genomes/mm10/hisat2index/mm10_spliceSites.txt"
  GENEMODEL="/scratch/tmp/a_toen03/Genomes/mm10/Gencode_M20/GRCm38_mm10_Ensembl_GeneModel.bed"
  fi

########################################################################################################################################

if [[ -e missing_tools.txt ]]; then rm missing_tools.txt; fi

function PathCheck {
  
  if [[ $(command -v $1 | wc -l) == 0 ]]; then 
    echo ${1} >> missing_tools.txt
    fi
  
}; export -f PathCheck

TOOLS=(hisat2 samtools read_distribution.py read_duplication.py parallel xargs)

for i in $(echo ${TOOLS[*]}); do
  PathCheck $i; done
  
if [[ -e missing_tools.txt ]] && [[ $(cat missing_tools.txt | wc -l | xargs) > 0 ]]; then
  echo '[ERROR] Tools missing in PATH -- see missing_tools.txt' && exit 1
  fi
 
########################################################################################################################################

function AlnQualControl {
  
  BASENAME=$1 
  
  ## Check if fastq files are present
  if [[ ! -e ${BASENAME}_1.fastq.gz ]] || [[ ! -e ${BASENAME}_2.fastq.gz ]]; then
    (>&2 paste -d " " <(echo '[ERROR]' 'Missing input files for' $BASENAME ' -- exit') <(date))
    exit 1
    fi
  
  (>&2 paste -d " " <(echo '[INFO]' 'HisatALN for' $1 'started on') <(date))
  
  hisat2 -p 12 \
    --rg-id ${BASENAME} \
    --known-splicesite-infile $3 \
    --summary-file ${1}_hisat2_report.log \
    -x $2 \
    -1 ${BASENAME}_1.fastq.gz \
    -2 ${BASENAME}_2.fastq.gz | \
  samblaster | \
  tee >(samtools flagstat - > ${BASENAME}_sorted.flagstat) | \
  samtools view -m 1G -@ 12 -b | \
  tee ${BASENAME}_sorted.bam | \
  samtools index - > ${BASENAME}_sorted.bam.bai
  
  ## QC:
  read_distribution.py -i ${BASENAME}_sorted.bam -r $4 > ${BASENAME}_sorted_readDistro.txt
  read_duplication.py -i ${BASENAME}_sorted.bam -o ${BASENAME} -q 20
  
  (>&2 paste -d " " <(echo '[INFO]' 'HisatALN for' $1 'ended on') <(date))
  
}; export -f AlnQualControl 

ls *_1.fastq.gz | awk -F "_1" '{print $1}' | sort -u | \
  parallel -j 4 "AlnQualControl {} $HISAT2_IDX $SPLICE_FILE $GENEMODEL 2> {}.log"
