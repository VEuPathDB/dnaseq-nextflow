#!/usr/bin/env bash

set -euo pipefail
cp genome.fa.gz hold.fa.gz
gunzip hold.fa.gz
for seq in \$(samtools view -H result_sorted.bam | grep '^@SQ' | cut -f 2); do echo \${seq#*SN:}; done > regions.txt
samtools faidx hold.fa \$(cat regions.txt) > genome_reordered.fa
samtools faidx genome_reordered.fa
