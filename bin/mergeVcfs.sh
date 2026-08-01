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

if (( $# == 0 )); then
  echo "ERROR: mergeVcfs.sh: no input VCFs -- mergeExperiments needs at least one per-strain VCF (check params.vcfFiles)" >&2
  exit 1
fi

norm=()
for vcf in "$@"; do
  n="${vcf%.vcf.gz}.norm.vcf.gz"
  bcftools annotate -x "INFO,FORMAT/GL,FORMAT/DPR" "$vcf" -Oz \
    | bcftools norm -m -any --multi-overlaps . -Oz -o "$n"
  bcftools index --tbi "$n"
  norm+=( "$n" )
done

if (( ${#norm[@]} == 1 )); then
  # Single strain: nothing to merge, and bcftools merge rejects one input. The
  # normalized file already IS the contract -- verified that bcftools merge does
  # not re-collapse rows that norm -m -any split within a single file, so this
  # is shape-identical to the n>=2 output minus the other samples' columns.
  mv "${norm[0]}" "$out"
else
  bcftools merge --merge both -O z -o "$out" "${norm[@]}"
fi
