#!/usr/bin/env bash

set -euo pipefail
cp genome.fa hold.fa
hisat2-build hold.fa genomeIndex
samtools faidx hold.fa
