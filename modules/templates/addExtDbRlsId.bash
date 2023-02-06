#!/usr/bin/env bash

set -euo pipefail
perl /usr/bin/addExtDbRlsId.pl \
    --variationFile $variationFile \
    --gusConfig $gusConfig \
    --extdb_spec $extDbSpec
    
