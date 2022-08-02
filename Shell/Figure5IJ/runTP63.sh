#!/bin/sh
#$ -l mem_free=50G
#$ -l h_vmem=50G
#$ -pe smp 8

peakFile="TE5_TP63vsInput_summits.chr.bed"
fileIndex="TE5_TP63_H3K27ac_WGBS"

export PATH=$PATH:/common/bermanblab/build/homer/bin

cd Data/Figure5IJ/
TE5_TP63=TE5_TP63vsInput.SubtractCPM.bw
TE5_H3K27Ac=TE5_H3K27AcvsInput.SubtractCPM.bw
TE5_WGBS=TE5.sorted.rmblackList.bw

/hpc/home/zhengy2/.local/bin/computeMatrix reference-point -S ${TE5_TP63} ${TE5_H3K27Ac} ${TE5_WGBS} \
   -R ${peakFile} -b 1500 -a 1500 --binSize 30 --sortRegions descend --numberOfProcessors 8 --outFileName ${fileIndex}.gz \
   --outFileNameMatrix ${fileIndex}.tab --outFileSortedRegions ${fileIndex}.sort.bed
