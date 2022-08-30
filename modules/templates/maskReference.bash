#!/usr/bin/env bash

set -euo pipefail
perl /usr/bin/maskGenome.pl \
  -p result.pileup \
  -f genome_reordered.fa.fai \
  -dc $params.minCoverage \
  -o masked.fa
fold -w 60 masked.fa > ${sampleName}.masked.fa
samtools faidx ${sampleName}.masked.fa
