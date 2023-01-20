#!/usr/bin/env bash

set -euo pipefail

writeStudyConfig.pl \
  --name $params.experimentName \
  --file $ploidyFile \
  --outputFile ${sampleName}_ploidyConfig.txt \
  --sourceIdType NASequence \
  --protocol Ploidy \
  --profileSetName "${sampleName} (CNV)"
