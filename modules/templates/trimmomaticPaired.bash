#!/usr/bin/env bash

set -euo pipefail
mateAEncoding=\$(<mateAEncoding)
java org.usadellab.trimmomatic.TrimmomaticPE \
  -trimlog trimLog.txt $sampleFile -\$mateAEncoding \
  -baseout sample ILLUMINACLIP:$adaptorsFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
