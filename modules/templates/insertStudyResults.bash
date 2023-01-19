#!/usr/bin/env bash

set -euo pipefail

ga ApiCommonWorkflow::Main::WorkflowSteps::InsertStudyResults \\
  --inputDir \'.\' \\
  --configFile \'$configFile\' \\
  --extDbRlsSpec \'$extDbSpec\' \\
  --studyName \'$sampleName\' \\
  --commit   

echo "DONE"
