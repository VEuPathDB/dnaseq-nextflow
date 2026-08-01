#!/usr/bin/env bash
# Characterization tests for bin/mergeVcfs.sh and bin/mergeCoverageBeds.sh.
#
# The contract under test: whatever the input arity, mergeVcfs.sh emits a
# merged.vcf.gz with freebayes INFO and FORMAT/GL,FORMAT/DPR stripped and
# multiallelic records split into biallelic rows; mergeCoverageBeds.sh emits a
# coverage.tsv with a chrom/start/end/<sample>... header over union intervals.
# n=1 must satisfy that contract identically to n>=2, and n=0 must fail loudly.
#
# No -e: several cases assert on non-zero exit codes.
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
MERGE_VCFS="$ROOT/bin/mergeVcfs.sh"
MERGE_BEDS="$ROOT/bin/mergeCoverageBeds.sh"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

fail() { echo "FAIL: $1"; exit 1; }

# Writes a single-sample VCF named $1.vcf.gz whose records are given as
# tab-joined "POS REF ALT INFO GT" lines on stdin.
mkvcf() {
  local name="$1" dir="$2"
  {
    printf '##fileformat=VCFv4.2\n'
    printf '##contig=<ID=chr1,length=100000>\n'
    printf '##INFO=<ID=AO,Number=A,Type=Integer,Description="alt obs">\n'
    printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">\n'
    # Number=. rather than freebayes' G/R: these fields get stripped anyway, and
    # a fixed arity would make bcftools complain about the multiallelic row for
    # reasons unrelated to what is under test.
    printf '##FORMAT=<ID=GL,Number=.,Type=Float,Description="likelihoods">\n'
    printf '##FORMAT=<ID=DPR,Number=.,Type=Integer,Description="depth per allele">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$name"
    while IFS=$'\t' read -r pos ref alt info gt; do
      printf 'chr1\t%s\t.\t%s\t%s\t.\t.\t%s\tGT:GL:DPR\t%s:-1,-2,-3:10,5\n' \
        "$pos" "$ref" "$alt" "$info" "$gt"
    done
  } > "$dir/$name.vcf"
  bgzip -f "$dir/$name.vcf"
}

# ---------------------------------------------------------------- mergeVcfs n=1

d="$TMP/vcf1"; mkdir -p "$d"
# A multiallelic SNP that norm -m -any must split, plus a plain SNP.
printf '100\tA\tG,T\tAO=5,6\t1/2\n200\tC\tT\tAO=9\t1/1\n' | mkvcf SoloStrain "$d"

( cd "$d" && "$MERGE_VCFS" merged.vcf.gz SoloStrain.vcf.gz ) \
  || fail "mergeVcfs.sh exited non-zero on a single input"

[ -s "$d/merged.vcf.gz" ] || fail "n=1 produced no merged.vcf.gz"

# bgzip-compressed and indexable (a plain-gzip or uncompressed file fails here)
tabix -f -p vcf "$d/merged.vcf.gz" 2>/dev/null \
  || fail "n=1 merged.vcf.gz is not bgzip-compressed / tabix-indexable"

rows="$(bcftools view -H "$d/merged.vcf.gz" | wc -l)"
[ "$rows" = "3" ] || fail "n=1 expected 3 rows (A>G, A>T, C>T), got $rows"

# every ALT is single-allele -- multiallelics were split
multi="$(bcftools view -H "$d/merged.vcf.gz" | cut -f5 | grep -c ',' || true)"
[ "$multi" = "0" ] || fail "n=1 left $multi multiallelic ALT field(s) unsplit"

# INFO stripped
badinfo="$(bcftools view -H "$d/merged.vcf.gz" | cut -f8 | grep -cv '^\.$' || true)"
[ "$badinfo" = "0" ] || fail "n=1 left INFO content on $badinfo row(s)"

# GL and DPR stripped from FORMAT
fmt="$(bcftools view -H "$d/merged.vcf.gz" | cut -f9 | sort -u | tr '\n' ',')"
case "$fmt" in
  *GL*) fail "n=1 left FORMAT/GL in place (FORMAT=$fmt)" ;;
esac
case "$fmt" in
  *DPR*) fail "n=1 left FORMAT/DPR in place (FORMAT=$fmt)" ;;
esac

# the sample column survives, named after the input
samples="$(bcftools query -l "$d/merged.vcf.gz" | tr '\n' ',')"
[ "$samples" = "SoloStrain," ] || fail "n=1 expected sample SoloStrain, got $samples"

echo "PASS: mergeVcfs.sh n=1"

# ---------------------------------------------------------------- mergeVcfs n=2

d="$TMP/vcf2"; mkdir -p "$d"
printf '100\tA\tG,T\tAO=5,6\t1/2\n' | mkvcf StrainA "$d"
printf '200\tC\tT\tAO=9\t1/1\n'     | mkvcf StrainB "$d"

( cd "$d" && "$MERGE_VCFS" merged.vcf.gz StrainA.vcf.gz StrainB.vcf.gz ) \
  || fail "mergeVcfs.sh exited non-zero on two inputs"

samples="$(bcftools query -l "$d/merged.vcf.gz" | sort | tr '\n' ',')"
[ "$samples" = "StrainA,StrainB," ] \
  || fail "n=2 expected samples StrainA,StrainB, got $samples"

rows="$(bcftools view -H "$d/merged.vcf.gz" | wc -l)"
[ "$rows" = "3" ] || fail "n=2 expected 3 rows (A>G, A>T, C>T), got $rows"

multi="$(bcftools view -H "$d/merged.vcf.gz" | cut -f5 | grep -c ',' || true)"
[ "$multi" = "0" ] || fail "n=2 left $multi multiallelic ALT field(s) unsplit"

echo "PASS: mergeVcfs.sh n=2"

# ---------------------------------------------------------------- mergeVcfs n=0

d="$TMP/vcf0"; mkdir -p "$d"
err="$( cd "$d" && "$MERGE_VCFS" merged.vcf.gz 2>&1 )"; rc=$?
[ "$rc" -ne 0 ] || fail "mergeVcfs.sh exited 0 with no input VCFs"
case "$err" in
  *"no input VCFs"*) : ;;
  *) fail "mergeVcfs.sh n=0 message did not name the missing input: $err" ;;
esac
[ ! -e "$d/merged.vcf.gz" ] || fail "mergeVcfs.sh n=0 created a merged.vcf.gz anyway"

echo "PASS: mergeVcfs.sh n=0"

# ---------------------------------------------------------- mergeCoverageBeds

mkbed() { printf '%b' "$2" | bgzip > "$1"; }

# n=1: header has exactly one sample column and data rows equal the input bed
d="$TMP/bed1"; mkdir -p "$d"
solo="chr1\t0\t100\t5.5\nchr1\t200\t300\t9\n"
mkbed "$d/SoloStrain_coverage.bed.gz" "$solo"

( cd "$d" && "$MERGE_BEDS" coverage.tsv SoloStrain_coverage.bed.gz ) \
  || fail "mergeCoverageBeds.sh exited non-zero on a single input"

hdr="$(head -1 "$d/coverage.tsv")"
[ "$hdr" = "$(printf 'chrom\tstart\tend\tSoloStrain')" ] \
  || fail "n=1 header wrong: $hdr"

cols="$(head -1 "$d/coverage.tsv" | awk -F'\t' '{print NF}')"
[ "$cols" = "4" ] || fail "n=1 expected 4 header columns (3 + 1 sample), got $cols"

printf '%b' "$solo" > "$TMP/expected.bed"
tail -n +2 "$d/coverage.tsv" > "$TMP/got.bed"
diff -u "$TMP/expected.bed" "$TMP/got.bed" \
  || fail "n=1 data rows differ from the single input bed"

echo "PASS: mergeCoverageBeds.sh n=1"

# n=2: unchanged unionbedg behavior, header widens by one column per sample
d="$TMP/bed2"; mkdir -p "$d"
mkbed "$d/StrainA_coverage.bed.gz" "chr1\t0\t100\t5\n"
mkbed "$d/StrainB_coverage.bed.gz" "chr1\t50\t150\t9\n"

( cd "$d" && "$MERGE_BEDS" coverage.tsv StrainA_coverage.bed.gz StrainB_coverage.bed.gz ) \
  || fail "mergeCoverageBeds.sh exited non-zero on two inputs"

hdr="$(head -1 "$d/coverage.tsv")"
[ "$hdr" = "$(printf 'chrom\tstart\tend\tStrainA\tStrainB')" ] \
  || fail "n=2 header wrong: $hdr"

# unionbedg splits the overlap into 0-50, 50-100, 100-150 and fills gaps with 0
rows="$(tail -n +2 "$d/coverage.tsv" | wc -l)"
[ "$rows" = "3" ] || fail "n=2 expected 3 union intervals, got $rows"

cols="$(tail -n +2 "$d/coverage.tsv" | awk -F'\t' 'NR==1{print NF}')"
[ "$cols" = "5" ] || fail "n=2 expected 5 data columns, got $cols"

echo "PASS: mergeCoverageBeds.sh n=2"

# n=0: fail loudly, write nothing
d="$TMP/bed0"; mkdir -p "$d"
err="$( cd "$d" && "$MERGE_BEDS" coverage.tsv 2>&1 )"; rc=$?
[ "$rc" -ne 0 ] || fail "mergeCoverageBeds.sh exited 0 with no input beds"
case "$err" in
  *"no input coverage BEDs"*) : ;;
  *) fail "mergeCoverageBeds.sh n=0 message did not name the missing input: $err" ;;
esac
[ ! -e "$d/coverage.tsv" ] || fail "mergeCoverageBeds.sh n=0 created a coverage.tsv anyway"

echo "PASS: mergeCoverageBeds.sh n=0"
