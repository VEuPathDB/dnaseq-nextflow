#!/usr/bin/env bash

set -euo pipefail
vcftools \
    --gzvcf varscan.concat.vcf.gz \
    --keep-only-indels \
    --out output \
    --recode
