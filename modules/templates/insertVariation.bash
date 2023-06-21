#!/usr/bin/env bash

set -euo pipefail

ga ApiCommonData::Load::Plugin::InsertVariant \\
  --extDbRlsSpec \'$extDbRlsSpec\' \\
  --variantFile \'$variantFile\' \\
  --commit   

echo "DONE"
