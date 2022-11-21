#!/usr/bin/env bash

set -euo pipefail
checkUniqueDeflines.pl \
    -i $params.fastaDir > check.txt
