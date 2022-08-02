#!/bin/sh

regionFile=$1

for wgbsfile in `ls Data/ESCA_rmblackList_cov5/*.all.cov5.sorted.rmblackList.bed`
do
   fileIndex=`basename ${wgbsfile}`
   fileIndex=${fileIndex%%.all.cov5.sorted.rmblackList.bed}
   bedtools intersect -a Data/MMSeekR_PMDs/${regionFile} -b ${wgbsfile} -wa -wb > Data/Figure7BD/PMD_methylation/${regionFile%%.bed}.${fileIndex}.bed
   java -jar meta/CalculateKeyValue.jar CalculateWindowMeanValue Data/Figure7BD/PMD_methylation/${regionFile%%.bed}.${fileIndex}.bed Data/Figure7BD/PMD_methylation/${regionFile%%.bed}.${fileIndex}.txt
   rm Data/Figure7BD/PMD_methylation/${regionFile%%.bed}.${fileIndex}.bed
done
