#!/usr/bin/env bash

set -euo pipefail

if [ $vcfCount -gt 1 ]; then

    bcftools merge -o result.vcf.gz -O z *.vcf.gz
    gunzip result.vcf.gz
    sed -i 's/\\%//g' result.vcf
    bgzip result.vcf
    
else

    cp *.vcf.gz result.vcf.gz
    gunzip result.vcf.gz
    sed -i 's/\\%//g' result.vcf
    bgzip result.vcf

fi
