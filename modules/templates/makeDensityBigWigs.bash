#!/usr/bin/env bash

set -euo pipefail
LC_COLLATE=C sort -k1,1 -k2,2n snpDensity.bed > sorted.snpDensity.bed
LC_COLLATE=C sort -k1,1 -k2,2n indelDensity.bed > sorted.indelDensity.bed
bedGraphToBigWig sorted.snpDensity.bed genome_reordered.fa.fai snpDensity.bw
bedGraphToBigWig sorted.indelDensity.bed genome_reordered.fa.fai indelDensity.bw
