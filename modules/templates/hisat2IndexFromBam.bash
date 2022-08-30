#!/usr/bin/env bash

set -euo pipefail
pwd 
touch genomeIndex.1.ht2
samtools faidx genome.fa
