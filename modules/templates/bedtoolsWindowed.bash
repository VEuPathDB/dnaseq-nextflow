#!/usr/bin/env bash

set -euo pipefail
bedtools coverage \
  -counts \
  -sorted \
  -g $genome \
  -a $windows \
  -b $resultSortedGatkBam > windowedCoverage.bed
