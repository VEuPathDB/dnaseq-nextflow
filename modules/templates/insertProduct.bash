#!/usr/bin/env bash

set -euo pipefail

ga ApiCommonData::Load::Plugin::InsertVariantProductSummary \\
  --variantProductFile \'$productFile\' \\
  --commit   

echo "DONE"
