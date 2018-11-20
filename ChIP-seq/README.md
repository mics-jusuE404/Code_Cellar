# ChIP-seq folder:

## ENCODE_download.sh
This script can be run to download all fastq files that are listed in a typical [ENCODE metadata file](https://www.encodeproject.org/metadata/type=Experiment&assay_term_name=ChIP-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_term_name=Karpas-422/metadata.tsv)
It will load all single-end and fwd-read paired-end files (to not mix single and paired, treat paired as single) and then rename the files,
which are like `ENSXXXXXX.fastq.gz` into `Biosample_Assay_FileAccession_ExperimentTarget_Biolrep<int>_TechRep<int>.fastq.gz`
