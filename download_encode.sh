#!/bin/bash

## Given a metadata.tsv file from encode in $(pwd), download the fastq files and rename with experimental target,
## the biological replicate number and accession number:

## Download paths from metadata:
grep 'fastq.gz' metadata.tsv | cut -f1 | while read p; do echo 'https://www.encodeproject.org/files/'${p}'/@@download/'${p}'.fastq.gz' >> download.txt;  done

## Download:
cat download.txt | parallel -j 6 "wget {}"

## Rename fastq in the form EXPERIMENTTARGET_ACCESSION_BIOLREP.fastqgz, skipping unlisted files:
for i in *.fastq.gz; do 
  if [[ $(grep $i metadata.tsv) == "" ]]; then continue; fi
  mv $i "$(grep $i metadata.tsv | cut -f13)_${i%.fastq.gz}_biolRep$(grep $i metadata.tsv | cut -f25).fastq.gz"
  done




##FOR PE:
#!/bin/bash

## Rename fastq in the form EXPERIMENTTARGET_ACCESSION_BIOLREP.fastqgz, skipping unlisted files:
for i in *.fastq.gz; do 
  
  ## Check if file is listed in metadata.tsv
  if [[ $(grep $i metadata.tsv) == "" ]]
    then
    echo '[INFO] skipping' "${i}" 'as not present in metadada.tsv' && continue
    fi

  ## Get new name:
  Biosample=$(awk -v VAR=${i} 'FS="\t" { if ($43 ~ VAR) print $0 }' metadata.tsv | cut -f7 | \
              awk '{gsub(" ", "");print}')
  
  Accession=$(awk -v VAR=${i} 'FS="\t" { if ($43 ~ VAR) print $0 }' metadata.tsv | cut -f4 | \
              awk '{gsub(" ", "");print}')
  
  if [[ $(awk -v VAR=${i} 'FS="\t" { if ($43 ~ VAR) print $0 }' metadata.tsv | cut -f36) == "1" && \
        $(awk -v VAR=${i} 'FS="\t" { if ($43 ~ VAR) print $0 }' metadata.tsv | cut -f35) == "paired-ended" ]]
    then
    N=1
    fi
  if [[ $(awk -v VAR=${i} 'FS="\t" { if ($43 ~ VAR) print $0 }' metadata.tsv | cut -f36) == "2" && \
        $(awk -v VAR=${i} 'FS="\t" { if ($43 ~ VAR) print $0 }' metadata.tsv | cut -f35) == "paired-ended" ]]
    then
    N=2
    fi  
  
  BiolRep=$(awk -v VAR=${i} 'FS="\t" { if ($43 ~ VAR) print $0 }' metadata.tsv | cut -f31)
  
  ## Construct new name:
  echo '[INFO]: Renaming' $i 'into' ${Biosample}"_"${Accession}"_BiolRep"${BiolRep}"_"${N}".fastq.gz"
  mv ${i} ${Biosample}"_"${Accession}"_BiolRep"${BiolRep}"_"${N}".fastq.gz"
  
  done
