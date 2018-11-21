# ChIP-seq folder:

## ENCODE_download.sh
This script can be run to download all fastq files that are listed in a typical [ENCODE metadata file](https://www.encodeproject.org/metadata/type=Experiment&assay_term_name=ChIP-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_term_name=Karpas-422/metadata.tsv).

It will loop through all metadata files in `$(pwd)` and loading all single-end and fwd-read paired-end files and then renaming them,
from like `ENSXXXXXX.fastq.gz` into `[Biosample-Experiment Type-Target-BiolRep-TechRep-File Accession]`.

## Fq2Bam_ChIPseqSE.sh
Script to align single-end `*fastq.gz`to hg38, resulting in unfiltered BAMs `*_raw.bam` and those without unmapped and duplicated reads `*_sorted.bam` plus indices, flagstats, CPM-normalized bigwigs and fastqc reports for each. Simply run `./Fq2Bam.sh` and it does the rest. Needs bwa, fastqc, samtools, sambamba and deeptools in PATH.

## Collapse_Replicates.sh
This one goes through the files and combined all technical replicates belonging to one biological replicate into one `fastq.gz`.
For example these two files 

<p>`FOO-LY1_ChIP-seq_H2AFZ_BiolRep1_TechRep1_ENCFF002AVF.fastq.gz`</p>
<p>`FOO_ChIP-seq_H2AFZ_BiolRep1_TechRep2_ENCFF002AYT.fastq.gz`</p>

will become
`FOO-LY1_ChIP-seq_H2AFZ_BiolRep1_combined.fastq.gz`
It also prints lists with all files used for combining and all those rejected due to read length.
