#!/bin/sh
export PATH=$PATH:/common/bermanblab/bin/
DMRFile=$1
for file in `ls Data/ESCA_rmblackList_cov5/*.all.cov5.sorted.rmblackList.bed`
do
  fileIndex=`basename ${file}`
  fileIndex=${fileIndex%%.all.cov5.sorted.rmblackList.bed}
  bedtools intersect -a ${file} -b ${DMRFile} -wa | sort | uniq | cut -f 4 | awk -v sample=${fileIndex} '{sum+=$1} END {print sample "\t" sum/NR}'
done
