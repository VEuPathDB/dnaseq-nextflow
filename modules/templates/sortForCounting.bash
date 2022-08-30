#!/usr/bin/env bash

set -euo pipefail
samtools sort -n result_sorted_gatk.bam > result_sortByName.bam
