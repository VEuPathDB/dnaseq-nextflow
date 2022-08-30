#!/usr/bin/env bash

set -euo pipefail
bcftools concat \
  -a \
  -o varscan.concat.vcf varscan.snps.vcf.gz varscan.indels.vcf.gz
