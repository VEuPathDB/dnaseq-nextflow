#!/usr/bin/env bash

set -euo pipefail
bedtools coverage \
  -a windows.bed \
  -b heterozygousSNPs.vcf \
  -sorted \
  -g genome.txt \
  -counts > heterozygousDensity.bed
