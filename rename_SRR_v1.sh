#!/bin/bash

## Rename every file in $(pwd) that starts with SRRxxxxx(...) given it is listed in accessions2name.txt

## Check if <accessions2name.txt> is present:
if [[ ! -f accessions2name.txt ]]; then 
  echo '[ERROR] accessions2name.txt is missing'
  paste <(echo '======>') <(echo 'accessions2name.txt must have $1=SRRxxxxx & $2=the new name')
  exit 1
  fi
  
##Check for duplicates in the new names:
if [[ $(cut -f2 accessions2name.txt | \
        awk NF | sort | uniq -c | \
        sed 's/^[ \t]*//;s/[ \t]*$//' | cut -d " " -f1 | sort -k1,1n | tail -n 1) > 1 ]]
  then 
  echo '[ERROR]: accessions2name.txt contains duplicate names which would overwrite some files -- exiting'
  exit 1
  fi

## Rename using gsub in awk:
for i in SRR*
  do
  
  ## Extract the leading SRR identifier:
  LEAD="$(echo "${i}" | awk -F "[._]" '{print $1}')"
  
  ## Check if listed in rename file, else exit:
  if [[ $(grep "${LEAD}" accessions2name.txt | wc -l) == 0 ]]; then 
    echo '"${i}" not listed in accessions2name.txt -- skip' && continue
    fi
  
  ## Grep the new name:
  GSUB=$(grep "${LEAD}" accessions2name.txt | cut -f2)
  
  ## Rename:
  echo '[INFO]: Renaming' "${i}" 'into' "$(echo "${i}" | awk -v lead=$LEAD -v replace=$GSUB '{gsub(lead,replace);print}')"
  
  mv ${i} $(echo "${i}" | awk -v lead=$LEAD -v replace=$GSUB '{gsub(lead,replace);print}')
  
  done >> renaming_report.log
  
  
