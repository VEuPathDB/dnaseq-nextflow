#!/usr/bin/env bash
 
set -euo pipefail

cp input.vcf.gz merge.vcf.gz
gunzip merge.vcf.gz
sed -i 's/\\%//g' merge.vcf
mv merge.vcf merged.vcf
gzip merged.vcf
