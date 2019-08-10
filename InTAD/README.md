## Scrips related to correlation analysis of ATAC-seq regions towards gene expression using InTAD from Bioconductor.

This folder contains scripts that were modified from the InTAD (Bioconductor) source to have better control over memory usage.
We noted that with default settings and a larger number of signal/gene pairs (~500.000 - 10^6) we ran out of memory in 
multicore-mode (even on the 1.5Tb node) while in single-core mode it takes a very long time (> overnight).
As we aimed to also run several permutations we decided to make some modifications to the original scripts to allow multicore processing
while still avoiding memory issues.

Doing so, `findCorrelations.R` was modified to perform a chunk-wise processing of the data. 
For this we make the `InTadSig` as usual but then process (default) 10000 signals ( =elements from `allIdGnY`)
at a time without q-value correction. Results are then printed to disk for every chunk.
Doing so we avoid memory overhead, clearing the memory (except the `InTadSig` object) after every chunk iteration.
It still must be run on the 1.5TB or 3TB memory node.

After this is finished one can load the disk-saved results and perform qvalue correction manually.
