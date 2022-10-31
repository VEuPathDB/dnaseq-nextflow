#!/usr/bin/env bash

set -euo pipefail
perl /usr/bin/calculatePloidy \
  --outputDir . \
  --fpkmFile $sampleFile 
  --ploidy $ploidy \
  --sampleName $sampleName \
  --taxonId  $taxonId \
  --geneFootprints footprints
perl /usr/bin/calculateGeneCNVs \
  --gusConfigFile gusConfig \
  --fpkmFile $sampleFile \
  --ploidy $ploidy \
  --outputDir . \
  --sampleName $sampleName \
  --taxonId $taxonId \
  --geneFootPrints footprints
