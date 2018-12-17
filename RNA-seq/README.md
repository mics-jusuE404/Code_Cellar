# Scripts to work with RNA-seq data

## MergeRename_ICGC.sh
Given a number of RNA-seq fastq files from ICGC (downloaded from the Collaboratory with `icgc-get`) that look like
`1f124e4fc6e1970b9c20048435b11c0a.RNA-PAIRED-tumor_4158726_RNA-SN758_0091_BC0YEEACXX_7_ATTCCT-R1.fastq.gz`
it `cat`s together all lane replicates per donor ID (here `4158726`) and per read direction (here `R1`).
If there is only one file per donor and R1/2, simply rename the file. In the end, output for the above file would be
`4158726_RNAseq_combined_1.fastq.gz`.

