#!/bin/sh

H3K36me2File=$1

cd /data/
computeMatrix scale-regions -S Figure4C/${H3K36me2File} --regionBodyLength 100000 -b 100000 -a 100000 -bs 1000 -R Figure4C/ESCA_regions.bed --skipZeros \
   --outFileSortedRegions Figure4C/${H3K36me2File%%.bw}.extend100k.sort.bed --sortRegions keep --numberOfProcessors 4 \
   --outFileName Figure4C/${H3K36me2File%%.bw}.extend100k.gz --outFileNameMatrix Figure4C/${H3K36me2File%%.bw}.extend100k.tab
