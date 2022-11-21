#!/usr/bin/env bash

set -euo pipefail

writeStudyConfig.pl \
  --outputFile ${sampleName}_geneCNVConfig.txt \
  --file geneCNVFile \
  --name "${sampleName}_geneCNV (CNV)" \
  --protocol geneCNV \
  --sourceIdType gene \
  --profileSetName "${sampleName} - GeneCNV"
