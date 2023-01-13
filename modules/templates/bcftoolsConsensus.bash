#!/usr/bin/env bash

set -euo pipefail
bcftools consensus \
  -I \
  -f $genomeMaskedFasta $varscanConcatVcfGz > cons.fa
