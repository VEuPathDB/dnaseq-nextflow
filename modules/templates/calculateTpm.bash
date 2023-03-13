#!/usr/bin/env bash

set -euo pipefail
makeTpmFromHtseqCountsCNV.pl \
  --geneFootprintFile $geneFootprintFile \
  --countFile $counts \
  --tpmFile out.tpm

