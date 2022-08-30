#!/usr/bin/env bash

set -euo pipefail
bedtools genomecov \
  -bg \
  -ibam result_sorted_gatk.bam \
  -g genome_reordered.fa.fai >coverage.bed
