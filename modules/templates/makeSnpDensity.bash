#!/usr/bin/env bash

set -euo pipefail
zcat varscan.snps.vcf.gz | bedtools coverage \
  -a windows.bed \
  -b stdin -sorted \
  -g genome.txt \
  -counts > snpDensity.bed
zcat varscan.indels.vcf.gz | bedtools coverage \
  -a windows.bed \
  -b stdin -sorted \
  -g genome.txt \
  -counts > indelDensity.bed
