#!/usr/bin/env bash

set -euo pipefail
perl /usr/bin/findValues.pl \
     -i output.recode.vcf \
     -s ${sampleName} \
     -o output.tsv
