#!/bin/bash

## Get cpm-normalized bigwig files with library size corrected by a TMM normalization factor.
## Script will treat all *_dedup.bam as belonging to the same experiment and will feed them to TMM.

usage(){
>&2 echo '
---------------------------------------------------------------------------------------------------------------------

---- Bigwig creation, normalization and optional averaging starting from BAM files.
---- Normalization method is TMM from edgeR.
---- USAGE: -h | --help         : Show this message
            -b                  : BAM files to process, e.g. "*.bam" or "A_dedup.bam B_dedup.bam" 
                                  Use -b not --bam as --bam has a yet-unsolved bug!
            -a | --atacseq      : Set to trigger ATAC-seq mode (=counting cutting sites)           [FALSE]
            -s | --pairedend    : Set to count fragments instead of reads                          [FALSE]
            -c | --onlytmm      : Only calculate TMM factors but do not create bigwigs             [FALSE]
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
OnlyCounts="FALSE"
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
   "--atacseq")   set -- "$@" "-a"   ;;   #
   "--pairedend") set -- "$@" "-s"   ;;   #
   "--extend")    set -- "$@" "-e"   ;;   #
   "--peaks")     set -- "$@" "-p"   ;;   #
   "--onlytmm")   set -- "$@" "-c"   ;;   #
   "--mean")      set -- "$@" "-m"   ;;   #
   "--wiggle")    set -- "$@" "-w"   ;;   #
   "--rscript")   set -- "$@" "-r"   ;;   #
   "--jobs")      set -- "$@" "-j"   ;;   #
   "--threads")   set -- "$@" "-t"   ;;   #
   *)             set -- "$@" "$arg" ;;   # 
 esac
done

#############################################################################################################################################

while getopts acsmb:e:p:w:e:j:t: OPT       
  do   
  case ${OPT} in                       
    b) BAMS="${OPTARG}"          ;;
    a) ATACseq="TRUE"            ;;   
    c) OnlyCounts="TRUE"         ;;
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

## Summary:
>&2 echo ''
>&2 echo '---------------------------------------------------------------------------------------------'
>&2 echo '[Info] Running with these parameters:'
>&2 echo '       -b          = ' ${BAMS}
>&2 echo '       --atacseq   = ' "${ATACseq}"
>&2 echo '       --pairedend = ' "${PairedEnd}"
>&2 echo '       --extend    = ' "${Extend}"
>&2 echo '       --onlytmm   = ' "${OnlyCounts}"
>&2 echo '       --peaks     = ' "${Peaks}"
>&2 echo '       --mean      = ' "${Mean}"
>&2 echo '       --wiggle    = ' "${WiggleTools}"
>&2 echo '       --rscript   = ' "${Rscript}"
>&2 echo '       --jobs      = ' "${Jobs}"
>&2 echo '       --threads   = ' "${Threads}"
>&2 echo '---------------------------------------------------------------------------------------------'
>&2 echo ''

export BAMS
export ATACseq
export OnlyCounts
export PairedEnd
export Extend
export Peaks
export Mean
export WiggleTools
export Rscript
export Jobs
export Threads

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
TOOLS=(${WiggleTools} ${Rscript} bamCoverage featureCounts bc wigToBigWig bg2bw bedtools mawk)

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
tmp.NF <- calcNormFactors(object = raw.counts, method = c("TMM"), doWeighting = FALSE)

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

fcBasic="featureCounts -a ${Peaks}.saf -F SAF -o ${Peaks}.saf.countmatrix -T ${Threads}"

## customize featureCounts based on assay and sequencing type:
if [[ $ATACseq == "TRUE" ]]; then
  eval "${fcBasic}" --read2pos 5 ${BAMS} 2> ${Peaks}.saf.countmatrix.log
elif [[ $PairedEnd == "TRUE" ]]; then
  eval "${fcBasic}" -p ${BAMS} 2> ${Peaks}.saf.countmatrix.log
else
  eval "${fcBasic}" ${BAMS} 2> ${Peaks}.saf.countmatrix.log
fi  

## feed to edgeR and obtain effetive normalization factors:
(>&2 paste -d " " <(echo '[Info]' 'Calculating TMM factors'))

## if a TMMfactors.txt already exist, move it to other name:
if [[ -e TMMfactors.txt ]]; then mv TMMfactors.txt TMMfactors_existing.txt; fi

if [[ -e "${Peaks}".saf.countmatrix ]] && [[ $(head "${Peaks}".saf.countmatrix | wc -l | xargs) > 0 ]]; then
  cat "${Peaks}".saf.countmatrix | ${Rscript} calculateTMM.R
fi

if [[ ! -e TMMfactors.txt ]] || (( $(head TMMfactors.txt | awk NF | wc -l | xargs) < 3 )); then
  (>&2 paste -d " " <(echo '[Error]' 'Size factor file is empty, absent or only one single entry'))
  exit 1
fi

if [[ ! -e tmp_chromSizes.txt ]]; then
    ls *dedup.bam | head -n 1 | while read p; do samtools idxstats $p | cut -f1,2 > tmp_chromSizes.txt; done
fi

if [[ ${OnlyCounts} == "TRUE" ]]; then
  (>&2 paste -d " " <(echo '[Info]' '--onlytmm is TRUE, therefore no bigwigs are created. Exiting.'))
  exit 0
  fi
  
#############################################################################################################################################

## Bigwig with deeptools using the normalization factors:
function BiggyWiggy {
  
  BAM=$1
  
  Out=$(echo ${BAM} | awk '{gsub("_dedup.bam", "_TMM.bigwig");print}')
  
  if [[ ! -e $BAM ]]; then
    (>&2 paste -d " " <(echo '[Error]' $BAM 'does not exist'))
    exit 1
  fi  
      
  ## extract size factor and get reciprocal:
  SF=$(bc <<< "scale=10; $(grep "${BAM}" TMMfactors.txt | cut -f2)^-1")
    
  if [[ $ATACseq == "FALSE" ]]; then
    ## if no index, do index =)  
    if [[ ! -e ${BAM}.bai ]]; then
      (>&2 paste -d " " <(echo '[Info]' 'Indexing' $BAM 'as no index found'))
      samtools index -@ ${Threads} ${BAM}
    fi
  
    ## the basic command for bamCoverage:
    Basic="bamCoverage --bam ${BAM} -o ${Out} -p ${Threads} -bs 1 --scaleFactor ${SF}"
    ## Extension of reads to fragmentss if paired-end data using TLEN:
    if [[ $PairedEnd == "TRUE" ]]; then
      eval "${Basic}" -e 2> ${BAM%.bam}_bamCoverage.log
    fi
  
    ## Extension of single end reads do user-defined frag. length
    if [[ $Extend != "FALSE" ]] && [[ $PairedEnd == "FALSE" ]]; then
      eval "${Basic}" -e ${Extend} 2> ${BAM%.bam}_bamCoverage.log
    fi
  
    ## single-end but no extention
    if [[ $Extend == "FALSE" ]]; then
      eval "${Basic}" 2> ${BAM%.bam}_bamCoverage.log
    fi
    
  fi
  
  ## ATAC-seq mode, reduce reads to cutting site shifted by +4/-5 using awk:
  if [[ $ATACseq == "TRUE" ]]; then
  
    bedtools bamtobed -i ${BAM} \
    | mawk 'OFS="\t"{if($6 == "+")print $1,$2+4,$2+5,".",".",$6}{if($6 == "-")print$1,$3-5,$3-4,".",".",$6}' \
    | sort -k1,1 -k2,2n -k3,3n -k6,6 -S5G --parallel=${Threads} \
    | bedtools genomecov -bg -i - -g tmp_chromSizes.txt \
    | mawk -v ASF=${SF} 'OFS="\t" {print $1,$2,$3,$4*ASF}' \
    | bg2bw -i /dev/stdin -c tmp_chromSizes.txt -o ${Out}
    
  fi  
  
}; export -f BiggyWiggy

(>&2 paste -d " " <(echo '[Info]' 'Creating bigwigs for each sample'))

ls *_dedup.bam | parallel -j "${Jobs}" "BiggyWiggy {}"

(>&2 paste -d " " <(echo '[Info]' 'Done creating bigwigs for each sample'))

#############################################################################################################################################

## If set, average replicates using wiggletools
function BiggyAverage {
  
  Basename=$1
  WiggleTools=$2
  
  Files=$(ls ${Basename}*rep*_TMM.bigwig)
  if (( $(echo $Files | xargs | wc -l) > 1 )); then 
    echo "[Info]: Only one rep for" $Files
    exit 0
  fi
  
  ## for a reason I do not understand the bigwigs that wigToBw spills out are 10x larger than the 
  ## exact same file produced by bedGraph2bw therefore use that tool:
  ${WiggleTools} mean ${Files} \
  | wigToBigWig /dev/stdin tmp_chromSizes.txt ${Basename}_mean_TMM_big.bigwig
  bigWigToBedGraph ${Basename}_mean_TMM_big.bigwig /dev/stdout \
  | bg2bw -i /dev/stdin -c tmp_chromSizes.txt -o ${Basename}_mean_TMM.bigwig && \
  rm ${Basename}_mean_TMM_big.bigwig

}; export -f BiggyAverage 

if [[ $Mean == "TRUE" ]]; then
  
  (>&2 paste -d " " <(echo '[Info]' 'Averaging bigwigs across replicates'))
  
  if (( $(ls *rep*_TMM.bigwig 2> /dev/null | wc -l | xargs) == 0 )); then 
    (>&2 paste -d " " <(echo '[Error]' 'No bigwigs found to be averaged!'))
    exit 1
  fi
  
  ls *rep*_TMM.bigwig | awk -F "_rep" '{print $1 | "sort -u"}' | parallel -j 8 "BiggyAverage {} ${WiggleTools}"
  
fi
