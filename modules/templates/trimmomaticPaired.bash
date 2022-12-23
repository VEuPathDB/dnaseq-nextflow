#!/usr/bin/env bash

set -euo pipefail
mateAEncoding=\$(<mateAEncoding)

if [ "$params.trimmomaticAdaptorsFile" = "NA" ]; then
    java org.usadellab.trimmomatic.TrimmomaticPE \
        -trimlog trimLog.txt $sampleFile -\$mateAEncoding \
        -baseout sample ILLUMINACLIP:/usr/local/bin/All_adaptors-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
else
    java org.usadellab.trimmomatic.TrimmomaticPE \
        -trimlog trimLog.txt $sampleFile -\$mateAEncoding \
        -baseout sample ILLUMINACLIP:$params.trimmomaticAdaptorsFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
fi




