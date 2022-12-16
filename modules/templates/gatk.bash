#!/usr/bin/env bash

set -euo pipefail
JARPATH="/usr/GenomeAnalysisTK.jar"
java -jar \$JARPATH \
  -I picard.bam \
  -R genome_reordered.fa \
  -T RealignerTargetCreator \
  -o forIndelRealigner.intervals 2>realaligner.err
java -jar \$JARPATH \
  -I picard.bam \
  -R genome_reordered.fa \
  -T IndelRealigner -targetIntervals forIndelRealigner.intervals \
  -o result_sorted_gatk.bam 2>indelRealigner.err
