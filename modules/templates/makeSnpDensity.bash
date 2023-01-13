#!/usr/bin/env bash

set -euo pipefail
zcat $varscanSnpsVcfGz | bedtools coverage \
  -a $windows \
  -b stdin -sorted \
  -g $genome \
  -counts > snpDensity.bed
zcat $varscanIndelsVcfGz | bedtools coverage \
  -a $windows \
  -b stdin -sorted \
  -g $genome \
  -counts > indelDensity.bed
