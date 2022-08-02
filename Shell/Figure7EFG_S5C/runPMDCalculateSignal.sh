#!/bin/sh
cat Shell/Figure7EFG_S5C/runPMDCalculateSignal.txt | xargs -iFile -P15 bash -c "File"
