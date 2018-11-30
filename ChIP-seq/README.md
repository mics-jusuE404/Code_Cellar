# ChIP-seq folder:

## ENCODE_download.sh
This script can be run to download all fastq files that are listed in a typical [ENCODE metadata file](https://www.encodeproject.org/metadata/type=Experiment&assay_term_name=ChIP-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_term_name=Karpas-422/metadata.tsv).

It will loop through all metadata files in `$(pwd)` and loads all single-end and fwd-read paired-end files followed by renaming them,
from like `ENSXXXXXX.fastq.gz` into `[Biosample-Experiment Type-Target-BiolRep-TechRep-File Accession]`.


## Collapse_Replicates.sh
This one goes through the files and combined all technical replicates belonging to one biological replicate into one `fastq.gz`.
For example these two files 

1. `FOO_ChIP-seq_H2AFZ_BiolRep1_TechRep1_ENCFF002AVF.fastq.gz`
2. `FOO_ChIP-seq_H2AFZ_BiolRep1_TechRep2_ENCFF002AYT.fastq.gz`

will become
`FOO_ChIP-seq_H2AFZ_BiolRep1_combined.fastq.gz`
It also prints lists with all files used for combining and all those rejected due to read length.

## Fq2Bam_ChIPseqSE.sh
Script to align single-end `*fastq.gz`to hg38, resulting in unfiltered BAMs `*_raw.bam` and those without unmapped and duplicated reads `*_sorted.bam` plus indices, flagstats, CPM-normalized bigwigs and fastqc reports for each. Simply run `./Fq2Bam.sh` and it does the rest. Needs bowtie, fastqc, samtools, sambamba in PATH. All reads are trimmed to 36bp to avoid mappability bias for the few ENCODE datasets that have 76bp reads.

## To do:
Add a flag option to the `Fq2Bam_ChIPseqSE.sh` script that allows to specify if bowtie is used for very short reads < 50bp or BWA mem for longer things and paired-end data.
