#!/usr/bin/env bash

set -euo pipefail
LC_COLLATE=C sort -k1,1 -k2,2n $coverageBed > sorted.coverage.bed
bedGraphToBigWig sorted.coverage.bed $genomeReorderedFastaIndex coverage.bw
