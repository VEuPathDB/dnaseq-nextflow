#!/usr/bin/env bash

set -euo pipefail
htseq-count \
  -f bam \
  -s no \
  -t CDS \
  -i gene_id \
  -a 0 $resultSortByNameBam $gtfFile > counts.txt
