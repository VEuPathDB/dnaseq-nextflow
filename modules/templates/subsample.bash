#!/usr/bin/env bash

set -euo pipefail
samtools index result_sorted.bam

# number of mapped reads is col 3 from idxstats
frac=\$( samtools idxstats result_sorted.bam | awk 'BEGIN {total=0} {total += \$3} END {frac=$params.maxNumberOfReads/total;  if (frac > 1) {print 1} else {print frac}}' )

# this will subsample fraction of mapped reads
if awk "BEGIN {exit !(\$frac >= 1)}"
 then
   ln -s result_sorted.bam result_sorted_ds.bam

else
   samtools view -b -s \$frac result_sorted.bam > result_sorted_ds.bam
fi
