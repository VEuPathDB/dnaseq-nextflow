#!/usr/bin/env bash

set -euo pipefail

ga ApiCommonData::Load::Plugin::InsertVariantAlleleSummary \\
  --variantAlleleFile \'$alleleFile\' \\
  --commit   

echo "DONE"
