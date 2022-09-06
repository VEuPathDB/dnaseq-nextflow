#!/usr/bin/env bash

set -euo pipefail
bcftools consensus \
  -I \
  -f genome_masked.fa varscan.concat.vcf.gz > cons.fa
