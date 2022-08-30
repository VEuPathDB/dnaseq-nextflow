#!/usr/bin/env bash

set -euo pipefail
samtools mpileup \
  -A \
  -f genome_reordered.fa \
  -B result_sorted_gatk.bam > result.pileup 2>pileup.err
