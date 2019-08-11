## Scripts to run InTAD analysis

Usage of `InTAD_RunAnalysis.R` assuming to work on a 3TB RAM node with 144 cores:

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
