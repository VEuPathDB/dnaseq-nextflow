#!/usr/bin/env bash

set -euo pipefail

ga ApiCommonData::Load::Plugin::InsertProducts \\
  --extDbRlsSpec \'$extDbRlsSpec\' \\
  --genomeExtDbRlsSpec \'$genomeExtDbRlsSpec\' \\
  --ProductFile \'$productFile\' \\
  --commit   

echo "DONE"
