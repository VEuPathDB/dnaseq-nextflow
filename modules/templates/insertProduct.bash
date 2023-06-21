#!/usr/bin/env bash

set -euo pipefail

ga ApiCommonData::Load::Plugin::InsertVariantProductSummary \\
  --variantProductFile \'$variantProductFile\' \\
  --commit   

echo "DONE"
