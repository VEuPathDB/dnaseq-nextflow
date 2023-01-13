#!/usr/bin/env bash

set -euo pipefail
mv $varscanConcatVcf ${sampleName}.concat.vcf
bgzip ${sampleName}.concat.vcf
tabix -fp vcf ${sampleName}.concat.vcf.gz
