#!/bin/bash

## Rename all files in a folder assuming that the first "." marks the beginning of the suffix:
## $1 is the SRRXXXXX without any suffix
## $2 is the new name, e.g. $1=SRRXXXXX $2=stemcell_rep1

## Check if renaming list is present:
if [[ ! -f accessions2name.txt ]]; then 
  echo '[ERROR] accessions2name.txt is missing'
  paste <(echo '======>') <(echo 'accessions2name.txt must have $1=SRR & $2=new name')
  exit 1
  fi

## Check if renaming list is without duplicate names:
if [[ $(cut -f2 accessions2name.txt | \
        awk NF | sort | uniq -c | \
        sed 's/^[ \t]*//;s/[ \t]*$//' | cut -d " " -f1 | sort -k1,1n | tail -n 1) > 1 ]]
  then 
  echo '[ERROR]: accessions2.txt contains duplicate names which would overwrite some files -- exiting'
  exit 1
  fi

## Rename:
for i in SRR*
  do
  
  ## Suffix = everything after the first "."
  SUFFIX=$(echo $i | cut -d "." -f2-)
  
  ## Prefix = everything before the first ".":
  if [[ $(echo $i | awk -F "_" '{print $1"\n"$2}' | awk NF | wc -l) > 1 ]]; then
    PREFIX=$(echo $i | awk -F "_" '{print $1}')
  else
    PREFIX=$(echo $i | awk -F "." '{print $1}')
  fi
  
  ## Middle = everything after the prefix and before the suffix:
  if [[ $(echo $i | awk -F "_" '{print $1"\n"$2}' | awk NF | wc -l) > 1 ]]; then
    MIDDLE=$(echo $i | awk -F "_" '{print $2}' | awk -F "." '{print "_"$1}')
    fi
  
  ## skip if file is not listed:
  if [[ $(grep "${PREFIX}" accessions2name.txt | wc -l) == 0 ]]; then 
    echo '"${i}" not in accessions2name.txt -- skip' && continue
    fi
   
  NEWNAME=$(paste -d "" <(echo $(grep "${PREFIX}" accessions2name.txt | cut -f2)) \
                        <(echo $"${MIDDLE}") \
                        <(echo '.') \
                        <(echo "${SUFFIX}") \
            )
  
  mv $i $NEWNAME
  
  done
