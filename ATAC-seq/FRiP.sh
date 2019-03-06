function FRiP {
  
  PEAKS=$1
  CUTSITES=$2
  
  if [[ -e $1 ]] && [[ -e $2 ]]; then
    TOTAL=$(pigz -c -d -p 4 $CUTSITES | wc -l)
    READ=$(bedtools intersect \
           -a $CUTSITES \
           -b <(bedtools slop -b 250 -g tmp_chromSizes.txt -i $PEAKS | sort -k1,1 -k2,2n)\
           -wa -sorted | wc -l)
     paste <(echo ${2%_cutsites.bed.gz}) <(bc <<< "scale=6;$READ/$TOTAL") 
     fi

}; export -f FRiP

ls *cutsites.bed.gz | \
  awk -F "_cutsites.bed.gz" '{print $1}' | \
  parallel "FRiP {}_summits.bed {}_cutsites.bed.gz" > FRiP_scores.txt
