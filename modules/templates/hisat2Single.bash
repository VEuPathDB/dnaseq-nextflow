#!/usr/bin/env bash

set -euo pipefail
mateAEncoding=\$(<mateAEncoding)
hisat2 --no-spliced-alignment \
  -k 1 \
  -p $params.hisat2Threads \
  -q --\$mateAEncoding \
  -x $hisat2_index \
  -U sample_1p \
    | samtools view -bS - \
    | samtools sort - > result_sorted.bam
