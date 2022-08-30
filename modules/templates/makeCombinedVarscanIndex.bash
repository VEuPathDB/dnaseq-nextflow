#!/usr/bin/env bash

set -euo pipefail
bgzip varscan.concat.vcf
tabix -fp vcf varscan.concat.vcf.gz
