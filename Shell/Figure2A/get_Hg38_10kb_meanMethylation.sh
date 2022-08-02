#!/bin/sh

#$ -l h_vmem=50G

#export PATH=$PATH:/common/bermanblab/bin/
#export PATH=$PATH:/common/bermanblab/java/jre1.8.0_31/bin
#export PATH=$PATH:/common/bermanblab/java/jdk1.8.0_31/bin

fileIndex=$1
bedtools intersect -a meta/hg38_10k.bed -b Data/ESCA_rmblackList_cov5/${fileIndex}.all.sorted.rmblackList.bed -wa -wb ｜ bedtools intersect -a - -b meta/Takai_Jones_from_Fei.hg38.bed -v > Data/Figure2A/hg38_10k.${fileIndex}.rmblackListCGI.bed
