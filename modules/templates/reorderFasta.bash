#!/usr/bin/env bash

set -euo pipefail
for seq in \$(samtools view -H $resultSortedBam | grep '^@SQ' | cut -f 2); do echo \${seq#*SN:}; done > regions.txt
samtools faidx $genomeFasta \$(cat regions.txt) > genome_reordered.fa
samtools faidx genome_reordered.fa
