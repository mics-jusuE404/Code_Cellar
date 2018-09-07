####################################################################################################################################
####################################################################################################################################

## Version last modified: 07.09.2018
## Alignment of paired-end data to hg38 on high-memory nodes. 6 jobs in parallel on largesmp nodes should be possible.
## Using seqtk to interleave mate files, BWA mem for alignment, samblaster for dup. marking, and sambamba for sorting and indexing.
## Usage ./script.sh BASENAME

## Run with a helper script that contains the #SBATCH parameters and :
## ls *_1.fastq.gz | awk -F "_" '{print $1}' | parallel -j 6 "./script.sh {}"

####################################################################################################################################
####################################################################################################################################

BASENAME=$1

## Truseq adapter:
ADAPTER1="AGATCGGAAGAGC"
ADAPTER2="AGATCGGAAGAGC"

####################################################################################################################################

BWA_IDX=/scratch/tmp/a_toen03/Genomes/hg38/bwa_index_noALT_withDecoy/hg38_noALT_withDecoy.fa

####################################################################################################################################
####################################################################################################################################

echo '#############################################################################################################################'
echo '[START]' $BASENAME 'on:' && date && echo ''

seqtk mergepe ${BASENAME}_1.fastq.gz ${BASENAME}_2.fastq.gz | \
  cutadapt -j 4 -a $ADAPTER1 -A $ADAPTER2 --interleaved -m 18 --max-n 0.1 - | \
  bwa mem -v 2 -R '@RG\tID:'${BASENAME}'_ID\tSM:'${BASENAME}'_SM\tPL:Illumina' -p -t 16 ${BWA_IDX} /dev/stdin | \
  samtools fixmate -m -@ 2 -O SAM - - | \
  samblaster --ignoreUnmated | \
  sambamba view -f bam -S -l 1 -t 4 -o /dev/stdout /dev/stdin | \
  sambamba sort -m 250G --tmpdir=./ -l 5 -t 24 -o ${BASENAME}_SortedRmdup.bam /dev/stdin  

echo '[END]' $BASENAME 'on:' && date
echo '#############################################################################################################################'

####################################################################################################################################
####################################################################################################################################
