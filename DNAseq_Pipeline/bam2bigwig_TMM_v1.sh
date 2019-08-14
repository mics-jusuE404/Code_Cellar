#!/bin/bash

## Get cpm-normalized bigwig files with library size corrected by a TMM normalization factor.
## Script will treat all *_dedup.bam as belonging to the same experiment and will feed them to TMM.

usage(){
echo '
---------------------------------------------------------------------------------------------------------------------

---- Bigwig creation, normalization and optional averaging starting from BAM files.
---- Normalization method is TMM from edgeR.
---- USAGE: -h | --help         : Show this message
            -b | --bams         : BAM files to process, e.g. "*.bam" or "A_dedup.bam B_dedup.bam" 
                                  Use -b not --bam as --bam has a yet-unsolved bug!
            -a | --atacseq      : Set to trigger ATAC-seq mode (=counting cutting sites)           [FALSE]
            -s | --pairedend    : Set to count fragments instead of reads                          [FALSE]
            -e | --extend       : Number of bp to extend reads to fragments.                       [0]                                        
            -p | --peaks        : Reference peak list to create the count matrix from              [no default]
            -m | --mean         : If set perform averaging of replicates                   
            -w | --wiggle       : Path to wiggletools executable if not the default one.
            -r | --rscript      : Path to Rscript executable if not the default one.
            -j | --jobs         : Number of parallel jobs                                          [4]
            -t | --threads      : Number of threads for bamCoverage                                [16]       
            
---------------------------------------------------------------------------------------------------------------------
'
}; if [[ -z "$1" ]] || [[ $1 == -h ]] || [[ $1 == --help ]]; then usage; exit; fi

#############################################################################################################################################

## Set defaults:

ATACseq="FALSE" 
PairedEnd="FALSE"
Extend="FALSE"
Mean="FALSE"
WiggleTools="$HOME/anaconda3_things/anaconda3/envs/wiggle/bin/wiggletools"
Rscript="$HOME/anaconda3_things/anaconda3/envs/R_env/bin/Rscript"
Jobs=4
Threads=16

#############################################################################################################################################

for arg in "$@"; do                         # for every arg in the commandline array ("$@")
 shift                                      # Shift by one so as to skip the script name
 case "$arg" in
   "--bams")      set -- "$@" "-b"   ;;   # 
   "--atacseq")   set -- "$@" "-a"   ;;   #
   "--pairedend") set -- "$@" "-s"   ;;   #
   "--extend")    set -- "$@" "-e"   ;;   #
   "--peaks")     set -- "$@" "-p"   ;;   #
   "--mean")      set -- "$@" "-m"   ;;   #
   "--wiggle")    set -- "$@" "-w"   ;;   #
   "--rscript")   set -- "$@" "-r"   ;;   #
   "--jobs")      set -- "$@" "-j"   ;;   #
   "--threads")   set -- "$@" "-t"  ;;  #
   *)             set -- "$@" "$arg" ;;   # 
 esac
done

#############################################################################################################################################

while getopts asmb:e:p:w:e:j:t: OPT       
  do   
  case ${OPT} in                       
    b) BAMS="${OPTARG}"          ;;
    a) ATACseq="TRUE"            ;;        
    s) PairedEnd="TRUE"          ;;
    e) Extend="${OPTARG}"        ;;
    p) Peaks="${OPTARG}"         ;;
    m) Mean="TRUE"               ;;
    w) WiggleTools="${OPTARG}"   ;;
    r) Rscript="${OPTARG}"       ;;
    j) Jobs="${OPTARG}"          ;;
    t) Threads="${OPTARG}"       ;;
  esac
done	
echo $BAMS
## Summary:
echo ''
echo '---------------------------------------------------------------------------------------------'
echo '[Info] Running with these parameters:'
echo '       --bams      = ' ${BAMS}
echo '       --atacseq   = ' "${ATACseq}"
echo '       --pairedend = ' "${PairedEnd}"
echo '       --extend    = ' "${Extend}"
echo '       --peaks     = ' "${Peaks}"
echo '       --mean      = ' "${Mean}"
echo '       --wiggle    = ' "${WiggleTools}"
echo '       --rscript   = ' "${Rscript}"
echo '       --jobs      = ' "${Jobs}"
echo '       --threads   = ' "${Threads}"
echo '---------------------------------------------------------------------------------------------'
echo ''

#############################################################################################################################################

## Check if required tools are in PATH and/or callable:

if [[ -e missing_tools.txt ]]; then rm missing_tools.txt; fi

## Function that checks if required tools are callable:
function PathCheck {
  
  if [[ $(command -v $1 | wc -l | xargs) == 0 ]]; then 
    echo ${1} >> missing_tools.txt
  fi
    
}; export -f PathCheck

## All tools:
TOOLS=(${WiggleTools} ${Rscript} bamCoverage featureCounts bc wigToBigWig)

## Loop through the tools and write to missing_tools.txt all those not in PATH / callable:
for i in $(echo ${TOOLS[*]}); do
  PathCheck $i; done
  
## If tools are missing write them to <missing_tools.txt>
if [[ -e missing_tools.txt ]] && [[ $(cat missing_tools.txt | wc -l | xargs) > 0 ]]; then
  echo '[ERROR] Tools missing in PATH, see missing_tools.txt' && exit 1; fi

#############################################################################################################################################

cat <<EOF > calculateTMM.R
## Rscript to calculate TMM size factors from a count matrix read from stdin,
## which is then used in the bigwig script to normalize bigwig counts.

## turn off scientific notation to avoid rounding errors:
options(scipen=99999)

## check for presence of required packages:
packageS <- c("data.table", "edgeR")
if (length(grep("FALSE", (packageS %in% rownames(installed.packages())))) > 0){
  stop("Package(s): ", packageS[which( packageS %in% rownames(installed.packages()) == "FALSE")], " are not installed!")
}

suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(data.table))

## read data:
raw.counts <- fread('cat /dev/stdin', skip = 1, header = T, data.table = F)

if (ncol(raw.counts) == 7) {
  write.table(data.frame(c("ONESAMPLE")), sep="\t", col.names = F, row.names = T,
              quote = F, file = "./sizeFactors.txt")
  stop("Only one sample present, so skipping calculation!")              
}

## expecting featureCounts format remove non-count columns
raw.counts <- raw.counts[,7:ncol(raw.counts)]

## normFactors:
tmp.NF <- calcNormFactors(object = raw.counts, method = c("TMM"))

## raw library size:
tmp.LS <- colSums(raw.counts)

## effective normalization factor => each count in the bigwig will be divided by this
if (file.exists("./TMMfactors.txt")) {
  Do.Append <- TRUE
} else Do.Append <- FALSE

suppressWarnings(
  write.table(x = data.frame(Sample    = names(tmp.LS),
                             TMMfactor = tmp.NF * ( tmp.LS / 1000000)),
              file = paste0("./TMMfactors.txt"),
              append = Do.Append, sep = "\t",
              col.names = !Do.Append, row.names = FALSE, quote = FALSE)
)
EOF

#############################################################################################################################################

## Make SAF format from peak list
(>&2 paste -d " " <(echo '[Info]' 'Writing SAF file'))

awk 'OFS="\t" {print $1"_"$2+1"_"$3, $1, $2+1, $3, "+"}' "${Peaks}" > "${Peaks}".saf

## Countmatrix:
(>&2 paste -d " " <(echo '[Info]' 'Creating count matrix'))

featureCounts -a "${Peaks}".saf -F SAF ${Read2Pos} ${paired} -o "${Peaks}".saf.countmatrix -T ${Threads} ${BAMS} 2> "${Peaks}".saf.countmatrix.log

## feed to edgeR and obtain effetive normalization factors:
(>&2 paste -d " " <(echo '[Info]' 'Calculating TMM factors'))

if [[ -e "${Peaks}".saf.countmatrix ]] && [[ $(head "${Peaks}".saf.countmatrix | wc -l | xargs) > 0 ]]; then
  cat "${Peaks}".saf.countmatrix | ${Rscript} calculateTMM.R
fi

if [[ ! -e TMMfactors.txt ]] || (( $(head TMMfactors.txt | awk NF | wc -l | xargs) < 3 )); then
  (>&2 paste -d " " <(echo '[Error]' 'Size factor file is empty, absent or only one single entry'))
  exit 1
fi  

#############################################################################################################################################

## Bigwig with deeptools using the normalization factors:
function BiggyWiggy {
  
  BAM=${1}
  Out=$(echo ${BAM} | awk '{gsub("_dedup.bam", "_TMM.bigwig");print}')
  Threads=${2}
  Extend=${3}
  
  if [[ $Extend == "FALSE" ]]
    then
    Extend=""
  else
    Extend="-e "${3}
  fi  
  
  if [[ ! -e $BAM ]]; then
    (>&2 paste -d " " <(echo '[Error]' $BAM 'does not exist'))
    exit 1
    fi
    
  if [[ ! -e ${BAM}.bai ]]; then
    (>&2 paste -d " " <(echo '[Info]' 'Indexing' $BAM 'as no index found'))
    samtools index -@ ${Threads} ${BAM}
  fi
  
  ## Get the scaling factor for the respective sample from TMMfactors.txt.
  ## ! Careful, the TMM itself is intended for division but --scaleFactor in bamCoverage multiplies,
  ##   so we have to do SF^(-1) first using bc (basic calculator):
  
  SF=$(bc <<< "scale=10; $(grep "${BAM}" TMMfactors.txt | cut -f2)^-1")
  
  ## bigwig (optional extension from reads to fragments)
  bamCoverage --bam ${BAM} -o ${Out} -p ${Threads} -bs 1 ${Extend} --scaleFactor ${SF} 2> ${BAM%.bam}_bamCoverage.log
  
}; export -f BiggyWiggy

(>&2 paste -d " " <(echo '[Info]' 'Creating bigwigs for each sample'))

ls *_dedup.bam | parallel -j "${Jobs}" "BiggyWiggy {} ${Threads} ${Extend}"

#############################################################################################################################################

## If set, average replicates using wiggletools
function BiggyAverage {
  
  Basename=$1
  WiggleTools=$2
  
  Files=$(ls ${Basename}*rep*_TMM.bigwig)
  
  ${WiggleTools} mean ${Files} | wigToBigWig stdin tmp_chromSizes.txt ${Basename}_mean_TMM.bigwig

}; export -f BiggyAverage 

(>&2 paste -d " " <(echo '[Info]' 'Averaging bigwigs across replicates'))

if [[ ! -e tmp_chromSizes.txt ]]; then
  ls *dedup.bam | head -n 1 | while read p; do samtools idxstats $p | cut -f1,2 > tmp_chromSizes.txt; done
fi

ls *rep*_TMM.bigwig | awk -F "_rep" '{print $1 | "sort -u"}' | parallel -j 8 "BiggyAverage {} ${WiggleTools}"