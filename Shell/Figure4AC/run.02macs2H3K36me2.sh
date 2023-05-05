#!/bin/sh
export PATH=$PATH:/hpc/home/zhengy2/.local/bin
export PATH=$PATH:/common/bermanblab/bin


ipFile=$1
controlFile=$2
outputName=$3

inputPath=/Volumes/Yuan_2T/WGBSpaper/Experiment/Process_files/01.Mapping
outputPath=/Volumes/Yuan_2T/WGBSpaper/Experiment/Process_files/02.MACS2

macs2 callpeak -t ${inputPath}/${ipFile}.sorted.mkdup.rmblacklist.bam -c ${inputPath}/${controlFile}.sorted.mkdup.rmblacklist.bam \
               -g 2.7e9 --broad -p 0.01 --extsize=146 --nomodel -n ${outputName} --outdir ${outputPath}

bamCompare -b1 ${inputPath}/${ipFile}.sorted.mkdup.rmblacklist.bam -b2 ${inputPath}/${controlFile}.sorted.mkdup.rmblacklist.bam --operation ratio \
           -o ${outputPath}/${outputName}.ratio_5k.bw --centerReads --binSize 5000 --scaleFactorsMethod None \
           --numberOfProcessors 2  --normalizeUsing CPM --ignoreDuplicates --extendReads 146