#!/bin/bash

## Joe Healey's extractor code, put into .bashrc:

function Extract {
  if [[ -f ${1} ]]
    then
    case ${1} in
      *.tar.bz2)   tar xvjf $1   ; echo "tar xvjf $1"  ;;
      *.tar.gz)    tar xvzf $1   ; echo "tar xvzf $1"  ;;
      *.bz2)       bunzip2 $1    ; echo "bunzip2 $1"   ;;
      *.rar)       unrar x $1    ; echo "unrar x $1"   ;;
      *.gz)        gunzip $1     ; echo "gunzip $1"    ;;
      *.tar)       tar xvf $1    ; echo "tar xvf $1"   ;;
      *.tbz2)      tar xvjf $1   ; echo "tar xvjf $1"  ;;
      *.tgz)       tar xvzf $1   ; echo "tar xvzf $1"  ;;
      *.zip)       unzip $1      ; echo "unzip $1"     ;;
      *.Z)         uncompress $1 ; echo "unzip $1"     ;;
      *.7z)        7z x $1       ; echo "7z x $1"      ;;
      *)           echo "dont know how to extract ${1}..."
    esac
  else
    echo "${1} is not a valid file!"
  fi
}
