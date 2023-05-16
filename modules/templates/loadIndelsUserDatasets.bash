#!/usr/bin/env bash

set -euo pipefail

ga ApiCommonData::Load::Plugin::InsertIndel \\
  --extDbRlsSpec \'$extDbRlsSpec\' \\
  --genomeExtDbRlsSpec \'$genomeExtDbRlsSpec\' \\
  --IndelFile \'$indel\' \\
  --commit   

cp $indel ${sampleName}_indel.tsv
  
echo "DONE"
