#!/usr/bin/env bash

set -euo pipefail

perl /usr/bin/writeStudyConfig.pl \
  --outputFile ${sampleName}_geneCNVConfig.txt \
  --file geneCNVFile \
  --name "${sampleName}_geneCNV (CNV)" \
  --protocol geneCNV \
  --sourceIdType gene \
  --profileSetName "${sampleName} - GeneCNV"
