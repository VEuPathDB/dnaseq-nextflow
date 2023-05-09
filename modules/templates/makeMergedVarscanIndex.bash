#!/usr/bin/env bash

set -euo pipefail
cp $resultVcfGz hold.vcf.gz
gunzip hold.vcf.gz
sed -i 's/\\%//g' hold.vcf
bgzip hold.vcf
mv hold.vcf.gz result.vcf.gz
tabix -fp vcf result.vcf.gz
