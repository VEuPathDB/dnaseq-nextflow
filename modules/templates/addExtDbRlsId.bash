#!/usr/bin/env bash

set -euo pipefail
addExtDbRlsId.pl \
    --variationFile $variationFile \
    --gusConfig $gusConfig \
    --extdb_spec $extDbSpec
    
