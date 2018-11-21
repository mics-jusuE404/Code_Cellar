# ChIP-seq folder:

## ENCODE_download.sh
This script can be run to download all fastq files that are listed in a typical [ENCODE metadata file](https://www.encodeproject.org/metadata/type=Experiment&assay_term_name=ChIP-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_term_name=Karpas-422/metadata.tsv).

It will loop through all metadata files in `$(pwd)` and loading all single-end and fwd-read paired-end files and then renaming them,
from like `ENSXXXXXX.fastq.gz` into `[Biosample-Experiment Type-Target-BiolRep-TechRep-File Accession]`.
Files belonging to one biological replicate can then easily be identified by removing the file accessions from the end like:
e.g. `ls *.gz | awk -F "ENC" '{print $1}'

## Fq2Bam_ChIPseqSE.sh
Script to align single-end `*fastq.gz`to hg38, resulting in unfiltered BAMs `*_raw.bam` and those without unmapped and duplicated reads `*_sorted.bam` plus indices, flagstats, CPM-normalized bigwigs and fastqc reports for each. Simply run `./Fq2Bam.sh` and it does the rest. Needs bwa, fastqc, samtools, sambamba and deeptools in PATH.
