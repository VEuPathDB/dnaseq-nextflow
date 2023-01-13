#!/usr/bin/env bash

set -euo pipefail
bedtools genomecov \
  -bg \
  -ibam $resultSortedGatkBam \
  -g $genomeReorderedFastaIndex >coverage.bed
