#!/bin/sh
type=$1
cd Data/
for file in `ls ESCA_rmblackList_cov5/*.all.cov5.sorted.rmblackList.bed`
do
  fileIndex=`basename ${file}`
  fileIndex=${fileIndex%%.all.cov5.sorted.rmblackList.bed}
  cut -f 4 ${file} | awk -v sample=${fileIndex} -v type=${type} '{sum+=$1} END {print sample "\t" type "\t" sum/NR}'
done
