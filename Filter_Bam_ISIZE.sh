#!/bin/bash

#### Extract read pairs from a paired-end BAM file that have a ISIZE (column 9) smaller than a used-supplied value.
#### Usage ----- ./Filter_Bam_ISIZE.sh in.bam 1000 ------ would get pairs with ISIZE smaller 1kb
#### Note: Script uses mawk (https://invisible-island.net/mawk/mawk.html) which is way faster than standard awk

sambamba view -f sam -h -t 4 $1 | \
  mawk -v LEN=$2 '{if ($9 <= LEN && $9 >= -(LEN) && $9 != 0 || $1 ~ /^@/) print $0}' | \
  sambamba view -f bam -h -o ${1%.bam}_isize${2}.bam -
