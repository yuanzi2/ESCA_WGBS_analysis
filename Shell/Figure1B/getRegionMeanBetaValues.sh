#!/bin/sh
regionFile=$1
echo ${regionFile}

cd Data/
for file in `ls ESCA_rmblackList_cov5/*.all.cov5.sorted.rmblackList.bed`
do
  fileIndex=`basename ${file}`
  fileIndex=${fileIndex%%.all.cov5.sorted.rmblackList.bed}
  bedtools intersect -a ${file} -b ${regionFile} -wa | uniq | cut -f 4 | awk -v sample=${fileIndex} '{sum+=$1} END {print sample "\t", sum/NR}'
done
