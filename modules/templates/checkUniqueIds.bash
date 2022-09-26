#!/usr/bin/env bash

set -euo pipefail
perl /usr/bin/checkUniqueDeflines.pl \
     -i $params.fastaDir > check.txt
