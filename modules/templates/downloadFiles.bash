#!/usr/bin/env bash

set -euo pipefail

perl /usr/local/bin/getFilesFromSra.pl --strain $strain --idList $idList

if [ -f "${strain}_2.fastq" ]; then
    export isPaired="true"
else
    export isPaired="false"
fi
