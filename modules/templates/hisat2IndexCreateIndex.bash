#!/usr/bin/env bash

set -euo pipefail
cp genome.fa.gz hold.fa.gz
gunzip hold.fa.gz
hisat2-build hold.fa genomeIndex
samtools faidx hold.fa
