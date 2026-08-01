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

if (( $# == 0 )); then
  echo "ERROR: mergeCoverageBeds.sh: no input coverage BEDs -- mergeExperiments needs at least one per-strain coverage BED (check params.coverageFiles)" >&2
  exit 1
fi

names=()
for f in "$@"; do
  names+=( "$(basename "$f" _coverage.bed.gz)" )
done

# IFS=$'\t' rather than the IFS='\t' this logic used while it lived in the
# Nextflow process body: there, Groovy turned \t into a real tab before bash saw
# it. In a standalone script bash reads '\t' literally as two characters, and
# "${names[*]}" would join the sample names with the first of them -- a
# backslash. Same output bytes as before, but no longer dependent on Groovy.
printf 'chrom\tstart\tend\t%s\n' "$(IFS=$'\t'; echo "${names[*]}")" > "$out"

if (( $# == 1 )); then
  # Single strain: a union of one file is the file, so there is nothing to do --
  # and unionbedg rejects one input anyway. modules/snp.nf emits these already
  # sorted, with sub-minCoverage gaps omitted rather than zero-filled, so
  # -filler 0 would have nothing to fill against either.
  zcat "$1" >> "$out"
else
  bedtools unionbedg -names "${names[@]}" -filler 0 -i "$@" >> "$out"
fi
