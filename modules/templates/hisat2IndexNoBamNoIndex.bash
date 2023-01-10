#!/usr/bin/env bash

set -euo pipefail
cp genome.fa hold.fa
TMP=$params.hisat2Index
FILES=\$TMP*
for f in \$FILES; do cp "\$f" "genomeIndex\${f#\$TMP}" ; done
samtools faidx hold.fa
