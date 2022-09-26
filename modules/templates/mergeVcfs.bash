#!/usr/bin/env bash

set -euo pipefail
bcftools merge \
      -o merged.vcf.gz \
      -O z *.vcf.gz
cp merged.vcf.gz toSnpEff.vcf.gz
gunzip toSnpEff.vcf.gz
