#!/usr/bin/env bash

set -euo pipefail

writeStudyConfig.pl \
  --outputFile ${sampleName}_geneCNVConfig.txt \
  --file $geneCNVFile \
  --name $params.experimentName \
  --protocol geneCNV \
  --sourceIdType gene \
  --profileSetName "${sampleName} - GeneCNV"
