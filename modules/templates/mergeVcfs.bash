#!/usr/bin/env bash

set -euo pipefail
for i in *.vcf.gz; do tabix \$i; done
bcftools merge \
      -o merged.vcf.gz \
      -O z *.vcf.gz
cp merged.vcf.gz toSnpEff.vcf.gz
gunzip toSnpEff.vcf.gz
