#!/usr/bin/env bash

set -euo pipefail
JARPATH="/usr/GenomeAnalysisTK.jar"
java -jar \$JARPATH \
  -I $picardBam \
  -R $genomeReorderedFasta \
  -T RealignerTargetCreator \
  -o forIndelRealigner.intervals 2>realaligner.err
java -jar \$JARPATH \
  -I $picardBam \
  -R $genomeReorderedFasta \
  -T IndelRealigner -targetIntervals forIndelRealigner.intervals \
  -o ${sampleName}.bam 2>indelRealigner.err
