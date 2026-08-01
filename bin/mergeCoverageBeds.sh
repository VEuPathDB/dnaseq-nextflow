#!/usr/bin/env bash
# Build the coverage.tsv contract consumed by bin/processSequenceVariations.jl:
# a "chrom start end <sample>..." header followed by union coverage intervals.
# Per-strain inputs are 4-column bgzipped BED (chrom start end mean_depth) from
# modules/snp.nf; the sample name is the basename minus _coverage.bed.gz.
#
# Usage: mergeCoverageBeds.sh OUT.tsv IN_coverage.bed.gz [IN_coverage.bed.gz ...]

set -euo pipefail

out="${1:-}"
if [[ -z "$out" ]]; then
  echo "ERROR: mergeCoverageBeds.sh: missing output path" >&2
  exit 1
fi
shift

names=()
for f in "$@"; do
  names+=( "$(basename "$f" _coverage.bed.gz)" )
done

printf 'chrom\tstart\tend\t%s\n' "$(IFS=$'\t'; echo "${names[*]}")" > "$out"

bedtools unionbedg -names "${names[@]}" -filler 0 -i "$@" >> "$out"
