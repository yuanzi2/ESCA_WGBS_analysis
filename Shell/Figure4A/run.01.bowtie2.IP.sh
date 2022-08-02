#!/bin/sh
export _JAVA_OPTIONS="-Xms2G -Xmx2G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=/Users/yueyuanzheng/temp/"
export PATH=$PATH:/Users/yueyuanzheng/bin/bin
thread=6

inputFileName=$1
inputFilePath=/Volumes/Yuan_2T/WGBSpaper/Experiment/Process_files/H3K36me2_data/
outputPath=/Volumes/Yuan_2T/WGBSpaper/Experiment/Process_files/01.Mapping/

genomeIndex="/Volumes/Yuan_2T/Reference/hg38/indexed_bowtie2/hg38"
blackList="/Volumes/Yuan_2T/Reference/hg38/hg38_blacklist.bed"

echo ${inputFileName}
bowtie2 --sensitive -p ${thread} -x ${genomeIndex} -1 ${inputFilePath}/${inputFileName}_1.fastq.gz -2 ${inputFilePath}/${inputFileName}_2.fastq.gz -S ${outputPath}/${inputFileName}.sam
samtools view -f 0x2 -q 10 -bS ${outputPath}/${inputFileName}.sam > ${outputPath}/${inputFileName}.bam
samtools sort -@ ${thread} -o ${outputPath}/${inputFileName}.sorted.bam ${outputPath}/${inputFileName}.bam
java -Xms2G -jar ~/Downloads/software/picard.jar MarkDuplicates I=${outputPath}/${inputFileName}.sorted.bam O=${outputPath}/${inputFileName}.sorted.mkdup.bam M=${outputPath}/${inputFileName}.sorted.mkdup.txt USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
bedtools intersect -v -a ${outputPath}/${inputFileName}.sorted.mkdup.bam -b ${blackList} > ${outputPath}/${inputFileName}.sorted.mkdup.rmblacklist.bam
samtools index ${outputPath}/${inputFileName}.sorted.mkdup.rmblacklist.bam
