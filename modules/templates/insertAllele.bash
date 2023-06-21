#!/usr/bin/env bash

set -euo pipefail

ga ApiCommonData::Load::Plugin::InsertVariantAlleleSummary \\
  --variantAlleleFile \'$variantAlleleFile\' \\
  --commit   

echo "DONE"
