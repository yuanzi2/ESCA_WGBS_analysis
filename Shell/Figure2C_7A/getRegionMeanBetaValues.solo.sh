#!/bin/sh
pmdFile=$1
echo ${pmdFile}

for file in `ls Data/ESCA_rmblackList_cov5/*.all.cov5.sorted.rmblackList.bed`
do
  fileIndex=`basename ${file}`
  fileIndex=${fileIndex%%.all.cov5.sorted.rmblackList.bed}
  bedtools intersect -a ${file} -b meta/solo_WCGW_hg38.bed -wa -f 1 | bedtools intersect -a - -b ${pmdFile} -wa | uniq | cut -f 4 | awk -v sample=${fileIndex} -v type="solo" '{sum+=$1} END {print sample "\t" type "\t" sum/NR}'
done
