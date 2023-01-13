#!/usr/bin/env bash

set -euo pipefail
LC_COLLATE=C sort -k1,1 -k2,2n $heterozygousDensityBed > sorted.heterozygousDensity.bed
bedGraphToBigWig sorted.heterozygousDensity.bed $genomeReorderedFastaIndex heterozygousDensity.bw
