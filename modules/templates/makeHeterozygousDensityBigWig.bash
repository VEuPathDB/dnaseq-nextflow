#!/usr/bin/env bash

set -euo pipefail
LC_COLLATE=C sort -k1,1 -k2,2n heterozygousDensity.bed > sorted.heterozygousDensity.bed
bedGraphToBigWig sorted.heterozygousDensity.bed genome_reordered.fa.fai heterozygousDensity.bw
