#!/bin/bash

#### Extract read pairs from a paired-end BAM file that have a ISIZE (column 9) smaller than a used-supplied value.
#### Usage ----- ./Filter_Bam_ISIZE.sh in.bam 1000 ------ would get pairs with ISIZE smaller 1kb
#### Note: Script uses mawk (https://invisible-island.net/mawk/mawk.html) which is way faster than standard awk

samtools view -h $1 | \
  mawk -v LEN=$2 '{if ($9 <= LEN && $9 >= -(LEN) && $9 != 0 || $1 ~ /^@/) print $0}' | \
  samtools view -bh -o ${1}_isize.bam -
