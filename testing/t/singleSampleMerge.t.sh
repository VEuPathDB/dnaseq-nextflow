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
