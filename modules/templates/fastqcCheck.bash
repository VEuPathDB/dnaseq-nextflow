#!/usr/bin/env bash

set -euo pipefail

if [ "$fromBam" = true ]; then

    touch mateAEncoding

else

    fastqc_check.pl $fastqc_output mateAEncoding

fi
