#!/usr/bin/env bash

set -euo pipefail
if [$params.picardJar -ne "NA"]; then
  JARPATH=$params.picardJar
else 
  JARPATH="/usr/picard/picard.jar"
fi
java -jar \$JARPATH AddOrReplaceReadGroups I=result_sorted_ds.bam O=picard.bam RGID=$sampleName RGSM=$sampleName RGLB=NA RGPL=NA RGPU=NA
java -jar \$JARPATH CreateSequenceDictionary R=genome_reordered.fa UR=genome_reordered.fa
java -jar \$JARPATH BuildBamIndex I=picard.bam
java -jar \$JARPATH CollectAlignmentSummaryMetrics R=genome_reordered.fa I=picard.bam O=summaryMetrics.txt

