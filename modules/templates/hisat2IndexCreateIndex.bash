#!/usr/bin/env bash

set -euo pipefail
hisat2-build genome.fa genomeIndex
samtools faidx genome.fa
