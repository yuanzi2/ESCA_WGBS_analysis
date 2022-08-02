#!/bin/sh
cd /data1/yueyuan/scATACseq/process
cat Shell/Figure7EFG_S5C/runDMRCalculateSignal.txt | xargs -iFile -P5 bash -c "File"
