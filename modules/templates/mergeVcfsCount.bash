#!/usr/bin/env bash

set -euo pipefail
echo $vcfCount
bcftools merge \
  -o result.vcf.gz \
  -O z *.vcf.gz
