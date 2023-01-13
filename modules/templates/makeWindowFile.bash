#!/usr/bin/env bash

set -euo pipefail
makeWindowedBed.pl \
  --samtoolsIndex $genomeReorderedFastaIndex \
  --winLen $winLen
