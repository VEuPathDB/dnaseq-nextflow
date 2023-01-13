#!/usr/bin/env bash

set -euo pipefail
bedtools coverage \
  -a $windows \
  -b $heterozygousSNPs \
  -sorted \
  -g $genome \
  -counts > heterozygousDensity.bed
