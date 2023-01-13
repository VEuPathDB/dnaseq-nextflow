#!/usr/bin/env bash

set -euo pipefail
JARPATH="/usr/picard/picard.jar"
java -jar \$JARPATH AddOrReplaceReadGroups I=$resultSortedDsBam O=picard.bam RGID=$sampleName RGSM=$sampleName RGLB=NA RGPL=NA RGPU=NA
java -jar \$JARPATH CreateSequenceDictionary R=$genomeReorderedFasta UR=$genomeReorderedFasta
java -jar \$JARPATH BuildBamIndex I=picard.bam
java -jar \$JARPATH CollectAlignmentSummaryMetrics R=$genomeReorderedFasta I=picard.bam O=summaryMetrics.txt

