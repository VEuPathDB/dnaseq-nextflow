#!/usr/bin/env bash

set -euo pipefail

perl /usr/local/bin/getFilesFromSra.pl --strain $strain --idList $idList
