#!/usr/bin/env bash

set -euo pipefail

ga ApiCommonData::Load::Plugin::InsertAlleles \\
  --extDbRlsSpec \'$extDbRlsSpec\' \\
  --genomeExtDbRlsSpec \'$genomeExtDbRlsSpec\' \\
  --AlleleFile \'$alleleFile\' \\
  --commit   

echo "DONE"
