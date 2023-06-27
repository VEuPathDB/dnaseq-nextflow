#!/usr/bin/env bash
 
set -euo pipefail

for i in *.vcf.gz; do cp $i $i.tmp.vcf.gz; gunzip $i.tmp.vcf.gz; bgzip $i.tmp.vcf; cp $i.tmp.vcf.gz $i; rm $i.tmp.vcf.gz; tabix $i; done
bcftools merge \
      -o merged.vcf.gz \
      -O z *.vcf.gz
cp merged.vcf.gz merge.vcf.gz
gunzip merge.vcf.gz
sed -i 's/\\%//g' merge.vcf
mv merge.vcf merged.vcf
