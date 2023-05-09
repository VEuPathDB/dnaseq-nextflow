#!/usr/bin/env bash

set -euo pipefail
gunzip $resultVcfGz
sed -i 's/\\%//g' *.vcf
bgzip *.vcf
tabix -fp vcf $resultVcfGz
