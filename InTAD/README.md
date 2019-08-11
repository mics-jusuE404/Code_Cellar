## Scripts to run InTAD analysis

Prediction of target genes from a set of ATAC-seq peaks based on Pearson correlation between
the peak and nearby TSS within the same TAD. Uses the InTAD framework (citation see below) from Bioconductor.
The `InTAD_RunAnalysis.R` script is a wrapper for this analysis, with options for shuffling samples 
and parallelization. For now it supports up to two different cell types for the analysis as this is what we had at max by now.

Usage of `InTAD_RunAnalysis.R` assuming to work on a 3TB RAM node:

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=144
#SBATCH --partition=largesmp
#SBATCH --mem=2600G
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --job-name=InTAD_permutations
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --output=rasmussen.log

$HOME/anaconda3_things/anaconda3/envs/InTAD/bin/Rscript \
InTAD_RunAnalysis.R \
--study Rasmussen \
--cores 36 \
--nshuff 100 \
--nshuffParallel 4 \
--rdata 190807_InTAD_PreparedData.Rdata
```
That would launch analysis using 100 shufflings (permutations) for this study,
with 36 cores per shuffle working on 4 shuffles in parallel. 
A run with ~600.000 analysis pairs it will very roughly need 1.5Tb, sometimes less, sometimes much more,
so do not run more than 4 jobs in parallel on the 3Tb node to avoid memory shortages.
If using the 1.5TB node go down to `--nshuffParallel 2` at max.
 
 
 ## Helper scripts
 `findCutoff_mclustBIC.R` is a script that estimates a cutoff for separating expressed from non-expressed genes/transcripts.
 Only needs a count matrix as input, by default will perform `log2+1`. 
 Script is modified from the `filterGeneExpr()` function of InTAD. 
 By default the `InTAD_RunAnalysis.R` will than only keep transcripts for which at least one group of cells has 
 a group-avagare above this threshold. 
 
 `makeUniqueCombinations.R` is a script that will take the colnames from the signals and expression objects to be used
 in the InTAD analysis, spilling out all possible combinations (currently for up to two different cell types).
 This is necessary as in our project we do not have "exact-match" data (so RNA-seq and ATAC-seq from the same donor)
 but from different littermate mice. Therefore we shuffle the data and run several permutations to ensure significant
 results are trushworthy and not a product of a lucky combination of two samples.


## Reference
Okonechnikov K, Erkek S, Korbel JO, Pfister SM, Chavez L (2019):
InTAD: chromosome conformation guided analysis of enhancer target genes.” BMC bioinformatics, 20, 60. doi: 10.1186/s12859-019-2655-2

https://bioconductor.org/packages/release/bioc/html/InTAD.html
