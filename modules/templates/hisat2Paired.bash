#!/usr/bin/env bash

set -euo pipefail
mateAEncoding=\$(<mateAEncoding)
hisat2 --no-spliced-alignment \
  -k 1 \
  -p $params.hisat2Threads \
  -q --\$mateAEncoding \
  -x $hisat2_index \
  -1 sample_1p \
  -2 sample_2p  \
    | samtools collate -@ $params.samtoolsThreads -o output.bam -
samtools fixmate -@ $params.samtoolsThreads -m output.bam fix.bam
samtools sort -@ $params.samtoolsThreads -o sort.bam fix.bam
samtools markdup -@ $params.samtoolsThreads -r sort.bam result_sorted.bam
