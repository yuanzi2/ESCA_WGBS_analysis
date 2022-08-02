#!/bin/sh

submitFile=$1
motifFile=$2
TF=$3

cd Data/Figure5IJ
awk '{print $1 "\t" $2-1500 "\t" $3+1500}' ${submitFile} > ${submitFile%%_summits.chr.bed}_summits.chr.3kb.bed
annotatePeaks.pl ${submitFile%%_summits.chr.bed}_summits.chr.3kb.bed hg38 -noann -m ${motifFile} -cpu 4 > ${submitFile%%_summits.chr.bed}_summits.chr.3kb.${TF}.txt
