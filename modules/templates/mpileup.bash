#!/usr/bin/env bash

set -euo pipefail
samtools mpileup \
  -A \
  -f $genomeReorderedFasta \
  -B $resultSortedGatkBam > result.pileup 2>pileup.err
