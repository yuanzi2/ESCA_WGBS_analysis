#!/bin/sh
#$ -l mem_free=30G
#$ -l h_vmem=50G
#$ -pe smp 1

fileIndex=$1
savePath=$2
export PATH=$PATH:/common/bermanblab/build/homer/bin

Eso26=Data/H3K27ac/EAC_cells/Eso26_H3K27AcvsInput.SubtractCPM.bw
Flo1=Data/H3K27ac/EAC_cells/Flo1_H3K27AcvsInput.SubtractCPM.bw
JH=Data/H3K27ac/EAC_cells/JH_H3K27AcvsInput.SubtractCPM.bw
OACp4C=Data/H3K27ac/EAC_cells/OACp4C_H3K27AcvsInput.SubtractCPM.bw
OE19=Data/H3K27ac/EAC_cells/OE19_H3K27AcvsInput.SubtractCPM.bw
OE33=Data/H3K27ac/EAC_cells/OE33_H3K27AcvsInput.SubtractCPM.bw
SKGT4=Data/H3K27ac/EAC_cells/SKGT4_H3K27AcvsInput.SubtractCPM.bw

KYSE70=Data/H3K27ac/ESCC_cells/KYSE70_H3K27AcvsInput.SubtractCPM.bw
KYSE140=Data/H3K27ac/ESCC_cells/KYSE140_H3K27AcvsInput.SubtractCPM.bw
KYSE150=Data/H3K27ac/ESCC_cells/KYSE150_H3K27AcvsInput.SubtractCPM.bw
KYSE180=Data/H3K27ac/ESCC_cells/KYSE180_H3K27AcvsInput.SubtractCPM.bw
KYSE200=Data/H3K27ac/ESCC_cells/KYSE200_H3K27AcvsInput.SubtractCPM.bw
TE5=Data/H3K27ac/ESCC_cells/TE5_H3K27AcvsInput.SubtractCPM.bw
TE7=Data/H3K27ac/ESCC_cells/TE7_H3K27AcvsInput.SubtractCPM.bw

/hpc/home/zhengy2/.local/bin/computeMatrix scale-regions -S  ${Eso26} ${Flo1} ${JH} ${OACp4C} ${OE19} ${OE33} ${SKGT4} \
   ${KYSE70} ${KYSE140} ${KYSE150} ${KYSE180} ${KYSE200} ${TE5} ${TE7} \
   -R Data/Figure5CD_S3GH/${fileIndex}.bed --skipZeros --numberOfProcessors 4 --outFileName ${savePath}/${fileIndex}.ESCACellsH3K27ac.SubtractCPM.gz \
   --outFileNameMatrix ${savePath}/${fileIndex}.ESCACellsH3K27ac.SubtractCPM.tab --outFileSortedRegions ${savePath}/${fileIndex}_sort.ESCACellsH3K27ac.SubtractCPM.bed
