#!/usr/bin/env bash

set -euo pipefail
addSampleToDefline.pl \
     -i $consFasta \
     -o unique_ids.fa \
     -s $sampleName
gzip unique_ids.fa
