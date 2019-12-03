# Scripts to work with RNA-seq data

## salmon_quantification_v1.sh
Shell script to run salmon (v.1.0.0 and later) on all fastq files in `$(pwd)`.
Options:
```sh
-h | --help      : Show this message
-i | --idx       : the transcriptome index
-m | --mode      : single or paired-end data                       [single,paired]
-n | --noLength  : turn off length correction if 3`-enriched
                   libraries are quantified                       
-t | --threads   : number of threads per sample                    [8]
-j | --jobs      : number of parallel jobs                         [8]
-l | --libtype   : library type                                    [A]
-f | --fastqc    : guess what this does                            [TRUE]
-u | --trim      : use cutadapt to trim adapters                   [FALSE]
