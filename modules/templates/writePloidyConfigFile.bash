#!/usr/bin/env bash

set -euo pipefail

perl /usr/bin/writeStudyConfig.pl \
  --name "${sampleName}_ploidy (CNV)" \
  --file ploidyFile \
  --outputFile ${sampleName}_ploidyConfig.txt \
  --sourceIdType NASequence \
  --protocol Ploidy \
  --profileSetName "${sampleName} (CNV)"
