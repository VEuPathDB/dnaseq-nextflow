#!/usr/bin/env bash

set -euo pipefail

if [ $vcfCount > 1 ]; then

    bcftools merge \
        -o result.vcf.gz \
        -O z *.vcf.gz
    
else

    cp *.vcf.gz result.vcf.gz

fi
