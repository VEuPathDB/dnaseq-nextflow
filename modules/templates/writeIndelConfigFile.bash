#!/usr/bin/env bash

set -euo pipefail

writeStudyConfig.pl \
  --name $params.experimentName \
  --file $indelFile \
  --outputFile ${sampleName}_indelConfig.txt \
  --sourceIdType NASequence \
  --protocol Indel \
  --profileSetName "${sampleName} (Indel)"
