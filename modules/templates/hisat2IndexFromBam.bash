#!/usr/bin/env bash

set -euo pipefail
pwd
cp genome.fa hold.fa
touch genomeIndex.1.ht2
samtools faidx hold.fa
