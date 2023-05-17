#!/usr/bin/env bash

set -euo pipefail

runGeneCNVAndPloidyQuery --organismAbbrev $organismAbbrev --geneSourceIdOrthologFile geneSourceIdOrthologFile --chrsForCalcsFile chrsForCalcsFile
