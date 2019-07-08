#!/bin/bash

## Wrapper to run meme::dream on a set of two fasta files (foreground-background) to scan for
## motif enrichments followed by matching the motifs against a motif database (here HOCOMOCO)


############################################################################################################

## Check all tools are in path:
TOOLS=(dreme tomtom)

## Function to stop script if tools not found:
function PathCheck {
  if [[ $(command -v $1 | wc -l) == 0 ]]; then 
    echo ${1} >> missing_tools.txt
    fi
}; export -f PathCheck

## Loop through tools:
for i in $(echo ${TOOLS[*]}); do
  PathCheck $i
  done
  
## If tools are missing write them to <missing_tools.txt>
if [[ -e missing_tools.txt ]] && [[ $(cat missing_tools.txt | wc -l | xargs) > 0 ]]; then
  echo '[ERROR] Tools missing in PATH -- see missing_tools.txt' && exit 1
fi

############################################################################################################

## DREAM/TOMTOM function:
function DREME {
  
  FOREGROUND=$1
  BACKGROUND=$2
  HOCOMOCO=$3
  
  if [[ $1 == $2 ]]; then 
    echo '[ERROR] Fore- and background are identical based on file name -- exiting'
    exit 1
    fi
  
  ## Run dreme
  dreme -o ./${1%.fa}_dreme_vs_${2%.fa} -p $1 -n $2 
  
  ## Now check motifs against hocomoco with tomtom:
  cd ./${1%.fa}_dreme_vs_${2%.fa}
  tomtom -o ./tomtom ./dreme.txt $HOCOMOCO
  mv ./tomtom/* . && rm -r ./tomtom
      
}; export -f DREME

./DREME $1 $2 $3

## Usage example: ./DREME goreground.fa background.fa /path/to/hocomoco/in/meme/format.meme
