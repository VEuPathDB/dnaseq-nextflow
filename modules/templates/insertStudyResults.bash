#!/usr/bin/env bash

set -euo pipefail

ga ApiCommonData::Load::Plugin::InsertStudyResults \\
  --inputDir \'.\' \\
  --configFile \'configFile\' \\
  --extDbSpec \'$extDbSpec\' \\
  --studyName \'$sampleName\' \\
  --commit   

echo "DONE"
