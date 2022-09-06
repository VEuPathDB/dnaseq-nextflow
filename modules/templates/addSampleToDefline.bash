#!/usr/bin/env bash

set -euo pipefail
perl /usr/bin/addSampleToDefline.pl \
     -i cons.fa \
     -o unique_ids.fa \
     -s $sampleName
gzip unique_ids.fa
