#!/usr/bin/env bash

set -euo pipefail
mkdir genome
mv $genesGtf genome/genes.gtf
mv $sequencesFa genome/sequences.fa
gzip genome/sequences.fa
cp /usr/bin/snpEff/snpEff.config .
java -jar /usr/bin/snpEff/snpEff.jar build -gtf22 -noCheckCds -noCheckProtein -v genome
java -Xmx4g -jar /usr/bin/snpEff/snpEff.jar genome $mergedVcf > merged.ann.vcf
