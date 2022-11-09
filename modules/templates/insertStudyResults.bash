#!/usr/bin/env bash

set -euo pipefail

ga ApiCommonData::Load::Plugin::InsertStudyResults \\
  --inputDir \'.\' \\
  --configFile \'configFile\' \\
  --extDbSpec \'$extDbSpec\' \\
  --studyName \'$studyName\' \\
  --platformExtDbSpec \'$platformExtDbSpec\' \\
  --commit   

echo "DONE"
