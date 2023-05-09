#!/usr/bin/env bash

set -euo pipefail
gunzip $resultVcfGz
sed -i 's/\\%//g' result.vcf
cp result.vcf holf.vcf
bgzip hold.vcf
rm result.vcf
mv hold.vcf result.vcf
tabix -fp vcf $resultVcfGz
