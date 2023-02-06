#!/usr/bin/env bash

set -euo pipefail
perl /usr/bin/addExtDbRlsIdToVariation.pl \
    --variationFile $variationFile \
    --gusConfig $gusConfig \
    --extdb_spec $extDbSpec
    
