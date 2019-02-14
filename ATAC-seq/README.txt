# ATACseq_lowlevel script

## Created files:
###*_raw.bam                      
=> sorted but completely unfiltered alignments
###*_rawbackup.bam
=> same as above but unsorted for later use with bedtools bamtobed
###tmp_chromSizes.txt
=> chromSizes file based on the BAM header
###*_dup.bam=>
sorted BAM excl. mapq < 20, non-prim. alignments, unmapped reads and non-primary chromosomes
