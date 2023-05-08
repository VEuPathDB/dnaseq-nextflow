#!/usr/bin/env bash
 
set -euo pipefail

for i in *.vcf.gz; do tabix \$i; done
bcftools merge \
      -o merged.vcf.gz \
      -O z *.vcf.gz
cp merged.vcf.gz merge.vcf.gz
gunzip merge.vcf.gz
sed -i 's/\%//g' merge.vcf
mv merge.vcf merged.vcf
