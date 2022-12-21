#!/usr/bin/env bash

set -euo pipefail
mateAEncoding=\$(<mateAEncoding)

if [ "$params.trimmomaticAdaptorsFile" = "NA" ]; then
    adaptorsFile="/usr/local/bin/All_adaptors-PE.fa"
else
    adaptorsFile="$params.trimmomaticAdaptorsFile"
fi

java org.usadellab.trimmomatic.TrimmomaticPE \
  -trimlog trimLog.txt $sampleFile -\$mateAEncoding \
  -baseout sample ILLUMINACLIP:$adaptorsFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
