#!/usr/bin/env bash

set -euo pipefail
bcftools concat \
  -a \
  -o varscan.concat.vcf $varscanSnpsVcfGz $varscanIndelsVcfGz
