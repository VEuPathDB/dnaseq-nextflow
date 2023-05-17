#!/usr/bin/env bash

set -euo pipefail
addExtDbRlsIdToVariation.pl \
    --variationFile $variationFile \
    --gusConfig $gusConfig \
    --extdb_spec "$extDbSpec"


