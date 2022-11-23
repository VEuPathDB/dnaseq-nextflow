#!/usr/bin/env bash

set -euo pipefail
mkdir genome
mv genes.gtf genome
mv sequences.fa.gz genome
cp /usr/bin/snpEff/snpEff.config .
java -jar /usr/bin/snpEff/snpEff.jar build -gtf22 -noCheckCds -noCheckProtein -v genome
java -Xmx4g -jar /usr/bin/snpEff/snpEff.jar genome merged.vcf > merged.ann.vcf
