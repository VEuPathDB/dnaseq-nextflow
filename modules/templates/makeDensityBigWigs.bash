#!/usr/bin/env bash

set -euo pipefail
LC_COLLATE=C sort -k1,1 -k2,2n $snpDensity > sorted.snpDensity.bed
LC_COLLATE=C sort -k1,1 -k2,2n $indelDensity > sorted.indelDensity.bed
bedGraphToBigWig sorted.snpDensity.bed $genomeReorderedFastaIndex snpDensity.bw
bedGraphToBigWig sorted.indelDensity.bed $genomeReorderedFastaIndex indelDensity.bw
