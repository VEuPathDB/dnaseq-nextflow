#!/usr/bin/env bash

set -euo pipefail
pwd
cp genome.fa.gz hold.fa.gz
gunzip hold.fa.gz
touch genomeIndex.1.ht2
samtools faidx hold.fa
