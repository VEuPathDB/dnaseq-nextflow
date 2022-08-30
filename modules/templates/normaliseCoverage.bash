#!/usr/bin/env bash

set -euo pipefail
# NOTE final processing requires querying the DB so can stay in ReFlow
normaliseCoverageCNV.pl \
  --bedFile windowedCoverage.bed \
  --summaryMetrics summaryMetrics.txt
