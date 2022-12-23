#!/usr/bin/env bash

set -euo pipefail
mv varscan.concat.vcf ${sampleName}.concat.vcf
bgzip ${sampleName}.concat.vcf
tabix -fp vcf ${sampleName}.concat.vcf.gz
