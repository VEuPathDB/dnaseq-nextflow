#!/usr/bin/env bash

set -euo pipefail
addSampleToDefline.pl \
     -i cons.fa \
     -o unique_ids.fa \
     -s $sampleName
gzip unique_ids.fa
