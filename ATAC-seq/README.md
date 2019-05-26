Lowlevel pipeline for ATAC-seq data with current settings optimized to run in parallel on four samples via `GNU parallel` on a (at least) 72-core node with > 80Gb RAM.

ATACseq_lowlevel_v12.sh accepts paired/single-end fastq or uBAM files depending on the MODE parameter in the header of the script.
