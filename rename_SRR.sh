#!/bin/bash

## Given a folder with fastq.gz files from the SRA (SRRXXXXX.fastq.gz or SRRXXXXX_1/2.fastq.gz), rename these files
## according to accessions2name.txt where 
## $1 is the SRRXXXXX without any _1, _2 or .fastq.gz extension and 
## $2 is the new name, e.g. $1=SRRXXXXX $2=stemcell_rep1

if [[ ! -f accessions2name.txt ]]; then 
  echo '[ERROR] accessions2name.txt is missing'
  paste <(echo '======>') <(echo 'accessions2name.txt must have $1=SRR & $2=new name')
  exit 1
  fi

for i in *.fastq.gz
  do
  
  SRR1=${i%.fastq*}
  SRR2=${SRR1%_*}
  
  ## skip if file is not listed:
  if [[ $(grep $SRR2 accessions2name.txt | wc -l) == 0 ]]; then echo '$i not in accessions2name.txt -- skip' && continue; fi
  
  ## Middle:
  if [[ $(echo $i | awk -F "_" '{print $1"\n"$2}' | awk NF | wc -l) > 1 ]]; then
    MIDDLE=$(echo $i | awk -F "_" '{print $2}' | awk -F ".fastq" '{print "_"$1}')
    fi
  
  NEWNAME=$(paste -d "" <(echo $(grep $SRR2 accessions2name.txt | cut -f2)) <(echo $MIDDLE) <(echo '.fastq.gz'))
  
  mv $i $NEWNAME
  
  done
