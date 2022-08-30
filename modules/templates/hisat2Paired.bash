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
    | samtools view -bS \
    | samtools sort \
    | samtools rmdup - result_sorted.bam
