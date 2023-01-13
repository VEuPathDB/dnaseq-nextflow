#!/usr/bin/env bash

set -euo pipefail
vcftools \
    --gzvcf $varscanConcatVcfGz \
    --keep-only-indels \
    --out output \
    --recode
