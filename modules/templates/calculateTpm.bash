#!/usr/bin/env bash

set -euo pipefail
makeTpmFromHtseqCountsCNV.pl \
  --geneFootprintFile geneFootprintFile \
  --countFile counts.txt \
  --tpmFile tpm.txt
#NOTE downstream processing from here requires querying DBs and occasional reloading - leave in ReFlow
