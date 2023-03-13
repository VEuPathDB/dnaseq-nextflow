#!/usr/bin/env bash

set -euo pipefail

calculatePloidy.pl \
    --outputDir . \
    --fpkmFile $sampleFile \
    --sampleName $sampleName \
    --taxonId  $taxonId \
    --geneFootprints $footprints \
    --ploidy $ploidy \
    --chrsForCalcsFile $chrsForCalc 

calculateGeneCNVs.pl \
    --fpkmFile $sampleFile \
    --ploidy $ploidy \
    --outputDir . \
    --sampleName $sampleName \
    --taxonId $taxonId \
    --geneFootPrints $footprints \
    --geneSourceIdOrthologFile $geneSourceIdOrtholog \
    --chrsForCalcsFile $chrsForCalc 


