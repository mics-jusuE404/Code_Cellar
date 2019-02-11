#!/bin/bash

## Quality control and alignment of PE RNA-seq data with hisat2, RSeQC, fastqc:

################################################################################################################################################

BASENAME=$1

## HISAT2 index, splice-site file and genome GTF:
HISAT2_IDX="/scratch/tmp/a_toen03/Genomes/mm10/hisat2_index/mm10"
SPLICE_FILE="/scratch/tmp/a_toen03/Genomes/mm10/hisat2index/mm10_spliceSites.txt"
GENEMODEL="/scratch/tmp/a_toen03/Genomes/mm10/Gencode_M20/GRCm38_mm10_Ensembl_GeneModel.bed"

function HisatALN {
  
  BASENAME=$1 
  
  ## Check if fastq files are present
  if [[ ! -e ${BASENAME}_1.fastq.gz ]] || [[ ! -e ${BASENAME}_2.fastq.gz ]]; then
    (>&2 paste -d " " <(echo '[ERROR]' 'Missing input files for' $BASENAME ' -- exit') <(date))
    exit 1
    fi
  
  (>&2 paste -d " " <(echo '[INFO]' 'HisatALN for' $1 'started on') <(date))
  
  ADAPTER1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
  ADAPTER2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
  cutadapt -j 12 -a $ADAPTER1 -A $ADAPTER2 -m 18 --max-n 0.1 \
    -o ${BASENAME}-trimmed-pair1.fastq.gz -p ${BASENAME}-trimmed-pair2.fastq.gz \
    ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz

  hisat2 -p 12 -X 1000 \
    --known-splicesite-infile $3 \
    --summary-file ${1}_hisat2_report.log \
    -x $2 \
    -1 ${BASENAME}-trimmed-pair1.fastq.gz \
    -2 ${BASENAME}-trimmed-pair2.fastq.gz | \
  samblaster | \
  sambamba view -S -f bam -l 5 -t 4 -o /dev/stdout /dev/stdin | \
  sambamba sort -m 4G --tmpdir=./ -l 6 -t 12 -o /dev/stdout /dev/stdin | \
  tee ${BASENAME}_sorted.bam | \
  tee >(samtools flagstat - > ${BASENAME}_sorted.flagstat) | \
  samtools index - > ${BASENAME}_sorted.bam.bai
  
  ## QC:
  read_distribution.py -i ${BASENAME}_sorted.bam -r $4 > ${BASENAME}_sorted_readDistro.txt
  read_duplication.py -i ${BASENAME}_sorted.bam -o ${BASENAME} -q 20
    
}; export -f HisatALN 

ls *.fastq.gz | parallel -j 72 "fastqc {}"

ls *_1.fastq.gz | awk -F "_1" '{print $1}' | sort -u | \
  parallel -j 4 "HisatALN {} $HISAT2_IDX $SPLICE_FILE $GENEMODEL 2>> {}.log"
  
  multiqc -o multiqc_all ./
