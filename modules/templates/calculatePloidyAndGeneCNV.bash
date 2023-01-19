#!/usr/bin/env bash

set -euo pipefail
calculatePloidy --outputDir . --fpkmFile $sampleFile --sampleName $sampleName --taxonId  $taxonId --geneFootprints footprints --ploidy $ploidy 
calculateGeneCNVs --gusConfigFile gusConfig --fpkmFile $sampleFile --ploidy $ploidy --outputDir . --sampleName $sampleName --taxonId $taxonId --geneFootPrints footprints
