#!/usr/bin/env bash

set -euo pipefail
JARPATH="/usr/local/VarScan.jar"
echo $sampleName >vcf_sample_name
java -jar \$JARPATH mpileup2snp result.pileup --vcf-sample-list vcf_sample_name --output-vcf 1 --p-value $params.varscanPValue --min-coverage $params.minCoverage --min-var-freq $params.varscanMinVarFreqSnp >varscan.snps.vcf  2>varscan_snps.err
java -jar \$JARPATH mpileup2indel result.pileup --vcf-sample-list vcf_sample_name --output-vcf 1 --p-value $params.varscanPValue --min-coverage $params.minCoverage --min-var-freq $params.varscanMinVarFreqCons >varscan.indels.vcf  2> varscan_indels.err
java -jar \$JARPATH mpileup2cns result.pileup --p-value $params.varscanPValue --min-coverage $params.minCoverage --min-var-freq $params.varscanMinVarFreqCons > varscan.cons 2>varscan_cons.err
bgzip varscan.snps.vcf
tabix -fp vcf varscan.snps.vcf.gz
bgzip varscan.indels.vcf
tabix -fp vcf varscan.indels.vcf.gz
bgzip varscan.cons
maskGenome.pl \
  -p result.pileup \
  -f genome_reordered.fa.fai \
  -dc $params.minCoverage \
  -o masked.fa
fold -w 60 masked.fa > genome_masked.fa
