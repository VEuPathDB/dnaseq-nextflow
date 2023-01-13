#!/usr/bin/env bash

set -euo pipefail
samtools sort -n $resultSortedGatkBam > result_sortByName.bam
