#!/usr/bin/env bash

set -euo pipefail
findValues.pl \
   -i $outputRecodeVcf \
   -s ${sampleName} \
   -o output.tsv
