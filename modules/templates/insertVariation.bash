#!/usr/bin/env bash

set -euo pipefail

ga ApiCommonData::Load::Plugin::InsertVariant \\
  --extDbRlsSpec \'$extDbRlsSpec\' \\
  --variantFile \'$variationFile\' \\
  --commit   

echo "DONE"
