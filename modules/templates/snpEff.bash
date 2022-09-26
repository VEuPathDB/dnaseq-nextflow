#!/usr/bin/env bash

set -euo pipefail
mkdir genome
mv genes.gtf.gz genome
mv sequences.fa.gz genome
cp /usr/bin/snpEff/snpEff.config .
if [ $params.databaseFileType = gtf ]; then
  java -jar /usr/bin/snpEff/snpEff.jar build -gtf22 -noCheckCds -noCheckProtein -v genome
elif [ $params.databaseFileType = gff ]; then
  java -jar /usr/bin/snpEff/snpEff.jar build -gff3 -noCheckCds -noCheckProtein -v genome
else
  echo "Params.databaseFileType is not gtf or gff"
fi
java -Xmx4g -jar /usr/bin/snpEff/snpEff.jar genome merged.vcf > merged.ann.vcf
