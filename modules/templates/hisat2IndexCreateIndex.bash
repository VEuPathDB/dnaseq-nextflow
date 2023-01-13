#!/usr/bin/env bash

set -euo pipefail
hisat2-build $genomeFasta genomeIndex
samtools faidx $genomeFasta
