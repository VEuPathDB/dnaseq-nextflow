#!/usr/bin/env bash

set -euo pipefail

if [ "$fromBam" = true ]; then

    touch sample_1P
    touch sample_2P
    
elif [ "$isPaired" = true ]; then
 
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
    
else

    touch sample_2P
    mateAEncoding=\$(<mateAEncoding)

    if [ "$params.trimmomaticAdaptorsFile" = "NA" ]; then
        java org.usadellab.trimmomatic.TrimmomaticSE \
            -trimlog trimLog.txt $sampleFile \
            -\$mateAEncoding sample_1P ILLUMINACLIP:/usr/local/bin/All_adaptors-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
    else
        java org.usadellab.trimmomatic.TrimmomaticSE \
            -trimlog trimLog.txt $sampleFile \
            -\$mateAEncoding sample_1P ILLUMINACLIP:$params.trimmomaticAdaptorsFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
    fi
    
fi
