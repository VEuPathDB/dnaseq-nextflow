#!/usr/bin/env bash

set -euo pipefail
JARPATH="/usr/local/VarScan.jar"
echo $sampleName >vcf_sample_name
java -jar \$JARPATH mpileup2snp $resultPileup --vcf-sample-list vcf_sample_name --output-vcf 1 --p-value $params.varscanPValue --min-coverage $params.minCoverage --min-var-freq $params.varscanMinVarFreqSnp >varscan.snps.vcf  2>varscan_snps.err
java -jar \$JARPATH mpileup2indel $resultPileup --vcf-sample-list vcf_sample_name --output-vcf 1 --p-value $params.varscanPValue --min-coverage $params.minCoverage --min-var-freq $params.varscanMinVarFreqCons >varscan.indels.vcf  2> varscan_indels.err
java -jar \$JARPATH mpileup2cns $resultPileup --p-value $params.varscanPValue --min-coverage $params.minCoverage --min-var-freq $params.varscanMinVarFreqCons > varscan.cons 2>varscan_cons.err
bgzip varscan.snps.vcf
tabix -fp vcf varscan.snps.vcf.gz
bgzip varscan.indels.vcf
tabix -fp vcf varscan.indels.vcf.gz
/usr/bin/perl /usr/bin/parseVarscanToCoverage.pl --file varscan.cons --percentCutoff 60 --coverageCutoff $params.minCoverage --outputFile ${sampleName}.coverage.txt
/usr/bin/perl /usr/bin/maskGenome.pl \
  -p $resultPileup \
  -f $genomeReorderedFastaIndex \
  -dc $params.minCoverage \
  -o masked.fa
fold -w 60 masked.fa > genome_masked.fa
