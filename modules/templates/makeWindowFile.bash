#!/usr/bin/env bash

set -euo pipefail
makeWindowedBed.pl \
  --samtoolsIndex genome_reordered.fa.fai \
  --winLen $winLen
