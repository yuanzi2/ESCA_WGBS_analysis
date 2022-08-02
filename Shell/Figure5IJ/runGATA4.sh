#!/bin/sh
#$ -l mem_free=50G
#$ -l h_vmem=50G
#$ -pe smp 8

peakFile="Eso26_Gata4flagvsInput_summits.chr.bed"
fileIndex="ESO26_GATA4_H3K27ac_WGBS"

export PATH=$PATH:/common/bermanblab/build/homer/bin

cd Data/Figure5IJ/
Eso26_GATA4=Eso26_Gata4flagvsInput.SubtractCPM.bw
Eso26_H3K27ac=ESO26_H3K27acvsInput.SubtractCPM.bw
Eso26_WGBS=ESO26.sorted.rmblackList.bw

/hpc/home/zhengy2/.local/bin/computeMatrix reference-point -S ${Eso26_GATA4} ${Eso26_H3K27ac} ${Eso26_WGBS} \
   -R ${peakFile}  -b 1500 -a 1500 --numberOfProcessors 8 --sortRegions descend --outFileName ${fileIndex}.gz --binSize 30 \
   --outFileNameMatrix ${fileIndex}.tab --outFileSortedRegions ${fileIndex}.sort.bed
