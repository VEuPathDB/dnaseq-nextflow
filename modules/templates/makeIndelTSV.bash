#!/usr/bin/env bash

set -euo pipefail
findValues.pl \
   -i output.recode.vcf \
   -s ${sampleName} \
   -o output.tsv
