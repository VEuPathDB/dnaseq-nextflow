#!/usr/bin/env bash

set -euo pipefail
LC_COLLATE=C sort -k1,1 -k2,2n coverage.bed > sorted.coverage.bed
bedGraphToBigWig sorted.coverage.bed genome_reordered.fa.fai coverage.bw
