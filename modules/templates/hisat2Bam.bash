#!/usr/bin/env bash

set -euo pipefail
samtools view -bS $sampleFile | samtools sort - > result_sorted.bam
