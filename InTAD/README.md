## Scrips related to correlation analysis of ATAC-seq regions towards gene expression using InTAD from Bioconductor.

This folder contains scripts that were modified from the InTAD source to have better control over memory usage.
We noted that with default settings and a larger number of signal/gene pairs (~500.000 - 10^6) we ran out of memory in 
multicore-mode (even on the 1.5Tb node) while in single-core mode it takes a very long time (> overnight).
As we aimed to also run several permutations we decided to make some modifications to the original scripts to allow multicore processing
while still avoiding memory issues.

Doing so, `findCorrelations.R` was modified to perform a chunk-wise processing of the data. 
For this we make the `InTadSig` as usual but then process 1000 signals ( =elements from `allIdGnY`)
at a time without q-value correction. Results are then printed to disk for every chunk.
Doing so we avoid memory overhead, clearing the memory from everything except the `InTadSig` after every chunk.
With 32 `mclapply`-workers this works now well on a standard 192GB RAM node.

After this is finished we then load the chunks that were saved to disk back into R using `data.table::fread`,
selecting only the relevant columns to perform q-value correction. 
