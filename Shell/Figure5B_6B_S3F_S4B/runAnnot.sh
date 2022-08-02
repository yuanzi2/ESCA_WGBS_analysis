#!/bin/sh
#$ -l mem_free=50G
#$ -l h_vmem=50G

export PATH=$PATH:~/bin/homer/bin

inputPath="Data/MaskUnionPMDs_DMRs/"
file=$1
annotatePeaks.pl ${inputPath}/${file} hg38 -cpu 4 > ${file%%.bed}.annot.gene.txt
