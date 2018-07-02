#!/bin/bash

##: Script to run lumpyexpress on a WGS tumor-normal pair that was aligned with BWA mem,
##: using the WGS_BWAmem.sh script.
##: "Raw" lumpy VCF will be piped to SVtyper to be genotyped.
##: Assumes that t/n-pairs are labelled ${BASENAME}-(t/n)_SortedRmdup.bam,
##: splitter and discordants are like ${BASENAME}-(t/n)_discordant.bam & ${BASENAME}-(t/n)_splitter.bam,
##: all BAMs are indexed as *.bam.bai
##: USAGE: ./lumpyexpress_SVtyper.sh Basename <single/paired>, e.g.:
## ./lumpyexpress_SVtyper.sh Sample1 paired

##: Written by Alexander Toenges (a.toenges@uni-muenster.de) version 05/2018
##########################################################################################################################################################################################
##########################################################################################################################################################################################

LUMPYEXP="$HOME/software/lumpy-sv-v0.2.13/bin/lumpyexpress"
SAMBAMBA="$HOME/software/sambamba"
SVTYPER="$HOME/software/svtyper-0.0.4/svtyper"
SNPEFF="java -jar $HOME/software/SnpEff/SnpSift.jar"

##########################################################################################################################################################################################
##########################################################################################################################################################################################

###########
BASENAME=$1
###########

## Check if $2 in <./lumpyexpress_SVtyper.sh Basename Layout> is either "single" for single-sample mode
## or "paired" for a tumor-normal pair:
if [[ $2 != "single" ]] && [[ $2 != "paired" ]]; then echo '$2 must either be <single> or <paired> -- exiting'; fi

##########################################################################################################################################################################################
##########################################################################################################################################################################################

##: PAIRED-MODE:
if [[ $2 == "paired" ]]
  then
  ##: Check if all files are present and indexed:
  ##: Full bam:
  if [[ ! -e ${BASENAME}-t_SortedRmdup.bam ]]; then echo '[ERROR]:' ${BASENAME}-t_SortedRmdup.bam 'does not seem to exist -- exiting' && exit 1; fi
  if [[ ! -e ${BASENAME}-n_SortedRmdup.bam ]]; then echo '[ERROR]:' ${BASENAME}-n_SortedRmdup.bam 'does not seem to exist -- exiting' && exit 1; fi

  ##: Splitter:
  if [[ ! -e ${BASENAME}-t_splitter.bam ]]; then echo '[ERROR]:' ${BASENAME}-t_splitter.bam 'does not seem to exist -- exiting' && exit 1; fi
  if [[ ! -e ${BASENAME}-n_splitter.bam ]]; then echo '[ERROR]:' ${BASENAME}-n_splitter.bam 'does not seem to exist -- exiting' && exit 1; fi

  ##: Discordants:
  if [[ ! -e ${BASENAME}-t_discordant.bam ]]; then echo '[ERROR]:' ${BASENAME}-t_discordant.bam 'does not seem to exist -- exiting' && exit 1; fi
  if [[ ! -e ${BASENAME}-n_discordant.bam ]]; then echo '[ERROR]:' ${BASENAME}-n_discordant.bam 'does not seem to exist -- exiting' && exit 1; fi

  ##: Indices:
  for i in ${BASENAME}*.bam
    do
      if [[ ! -e ${i}.bai ]]; then echo $i '[MAIN]:' 'not indexed -- indexing now:' &&  $SAMBAMBA index -t 16 $i; fi
    done

  ##: Lumpyexpress default mode with probability curve output for use with svtyper
  echo '[MAIN]: Running lumpyexpress/SVtyper on' $BASENAME 'tumor/normal pair'
  $LUMPYEXP -B ${BASENAME}-t_SortedRmdup.bam,${BASENAME}-n_SortedRmdup.bam -D ${BASENAME}-t_discordant.bam,${BASENAME}-n_discordant.bam -S ${BASENAME}-t_splitter.bam,${BASENAME}-n_splitter.bam -P -o ${BASENAME}_SVraw.vcf
  $SVTYPER -n 10000000 --bam ${BASENAME}-t_SortedRmdup.bam,${BASENAME}-n_SortedRmdup.bam > ${BASENAME}_SVgt.vcf
  
  ## Get somatic variants by selecting entries that have evidence in tumor but not in normal:
  $SNPEFF --file ${BASENAME}_SVgt.vcf "( isVariant( GEN[0] ) && ( GEN[1].AO == 0 ) )" | \
    mawk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > ${BASENAME}_SVgt_somatic.vcf
fi

##########################################################################################################################################################################################
##########################################################################################################################################################################################

if [[ $2 == "single" ]]
  then
  ## SINGLE-MODE:
  ##: Check if all files are present and indexed:
  ##: Full bam:
  if [[ ! -e ${BASENAME}_SortedRmdup.bam ]]; then echo '[ERROR]:' ${BASENAME}_SortedRmdup.bam 'does not seem to exist -- exiting' && exit 1; fi

  ##: Splitter:
  if [[ ! -e ${BASENAME}_splitter.bam ]]; then echo '[ERROR]:' ${BASENAME}_splitter.bam 'does not seem to exist -- exiting' && exit 1; fi

  ##: Discordants:
  if [[ ! -e ${BASENAME}_discordant.bam ]]; then echo '[ERROR]:' ${BASENAME}_discordant.bam 'does not seem to exist -- exiting' && exit 1; fi

  ##: Indices:
  for i in ${BASENAME}*.bam
    do
      if [[ ! -e ${i}.bai ]]; then echo $i '[MAIN]:' 'not indexed -- indexing now:' &&  $SAMBAMBA index -t 16 $i; fi
    done

  ##: Lumpyexpress default mode with probability curve output for use with svtyper
  echo '[MAIN]: Running lumpyexpress/SVtyper on' $BASENAME 'in single mode'
  $LUMPYEXP -B ${BASENAME}_SortedRmdup.bam -D ${BASENAME}_discordant.bam -S ${BASENAME}_splitter.bam -P -o ${BASENAME}_SVraw.vcf
  $SVTYPER -n 10000000 --bam ${BASENAME}_SortedRmdup.bam > ${BASENAME}_SVgt.vcf
fi    

##########################################################################################################################################################################################
##########################################################################################################################################################################################
