#!/usr/bin/env bash
# Combine per-strain VCFs into the merged.vcf.gz contract consumed by
# bin/processSequenceVariations.jl: freebayes INFO and FORMAT/GL,FORMAT/DPR
# stripped, and multiallelic records split into biallelic rows.
#
# Usage: mergeVcfs.sh OUT.vcf.gz IN.vcf.gz [IN.vcf.gz ...]

set -euo pipefail

out="${1:-}"
if [[ -z "$out" ]]; then
  echo "ERROR: mergeVcfs.sh: missing output path" >&2
  exit 1
fi
shift

norm=()
for vcf in "$@"; do
  n="${vcf%.vcf.gz}.norm.vcf.gz"
  bcftools annotate -x "INFO,FORMAT/GL,FORMAT/DPR" "$vcf" -Oz \
    | bcftools norm -m -any --multi-overlaps . -Oz -o "$n"
  bcftools index --tbi "$n"
  norm+=( "$n" )
done

bcftools merge --merge both -O z -o "$out" "${norm[@]}"
