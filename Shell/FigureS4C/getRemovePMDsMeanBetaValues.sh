#!/bin/sh
export PATH=$PATH:/common/bermanblab/bin/
type=$1
cd Data/
for file in `ls ESCA_cov3_maskUnionPMDs/*.all.sorted.rmblackList.rmPMDs.bed`
do
  fileIndex=`basename ${file}`
  fileIndex=${fileIndex%%.all.sorted.rmblackList.rmPMDs.bed}
  cut -f 4 ${file} | awk -v sample=${fileIndex} -v type=${type} '{sum+=$1} END {print sample "\t" type "\t" sum/NR}'
done
