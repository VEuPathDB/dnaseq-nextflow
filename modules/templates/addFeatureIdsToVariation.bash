#!/usr/bin/env bash

set -euo pipefail
addFeatureIdsToVariation.pl \
     --variationFile $variationFile \
     --gusConfig $gusConfig 

