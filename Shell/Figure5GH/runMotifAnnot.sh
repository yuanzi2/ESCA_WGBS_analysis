#!/bin/sh
#$ -pe smp 4
#$ -l mem_free=50G
#$ -l h_vmem=50G

export PATH=$PATH:/common/bermanblab/build/homer/bin
motifFile=CodeAndData/meta/HOCOMOCOv11_core_HUMAN_mono_homer_format_0.0001.motif
file=$1
saveFile=$2
annotatePeaks.pl ${file} hg38 -noann -m ${motifFile} -cpu 4 > ${saveFile}
