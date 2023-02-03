#!/usr/bin/env bash

set -euo pipefail

ga ApiCommonData::Load::Plugin::InsertSnps \\
  --extDbRlsSpec \'$extDbRlsSpec\' \\
  --genomeExtDbRlsSpec \'$genomeExtDbRlsSpec\' \\
  --SnpFile \'$snpFile\' \\
  --commit   

echo "DONE"
