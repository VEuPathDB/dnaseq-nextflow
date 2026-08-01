# Single-Sample mergeExperiments Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `mergeExperiments` work correctly when a single experiment contributes a single sample, without changing behavior for two or more samples.

**Architecture:** The shell bodies of `mergeVcfs` and `mergeCoverageBeds` move into `bin/mergeVcfs.sh` and `bin/mergeCoverageBeds.sh`, following the existing `bin/normalizePairedReadNames.sh` + `testing/t/normalizePairedReadNames.t.sh` pattern (Nextflow puts `bin/` on `PATH` automatically, so process bodies become one-line calls and the logic becomes testable from bash). Each script handles 1..N inputs internally and errors on 0. The workflow-level `single`/`multiple` branch is deleted, so the DAG no longer knows or cares about arity.

**Tech Stack:** Nextflow DSL2, bash, bcftools, bedtools, tabix/bgzip, Python 3 (comparator), all inside `veupathdb/dnaseqanalysis:latest`.

**Spec:** `docs/superpowers/specs/2026-08-01-single-sample-mergeexperiments-design.md`

**Branch:** `single-sample-merge` (already created; the spec commit is on it).

---

## Background the engineer needs

**What `mergeVcfs` produces is a contract, not just "a merge."** Downstream,
`bin/processSequenceVariations.jl` consumes `merged.vcf.gz` expecting freebayes
`INFO` and `FORMAT/GL`/`FORMAT/DPR` stripped and multiallelic records split into
biallelic rows. Both come from the per-file normalization loop, not from
`bcftools merge`. The current single-sample bypass in
`workflows/mergeExperiments.nf:29-31` skips that loop entirely — that is the bug.

**Two facts were verified in the container; do not re-litigate them:**

1. `bcftools merge` exits 1 on one input; `bedtools unionbedg` exits 1 on one
   input ("Only a single BedGraph file was specified").
2. `bcftools merge` does **not** re-collapse rows that `bcftools norm -m -any`
   split within a single file. A multiallelic `A -> G,T` split into two biallelic
   rows stays two rows. So for n=1, using the normalized file directly is
   shape-identical to what merge would have produced.

**Per-strain coverage BEDs are 4-column** (`chrom start end mean_depth`),
bgzipped, produced by `modules/snp.nf`. For a single file, `bedtools unionbedg`
output is byte-identical to the file's own contents, so `zcat` is exact.

**How to run anything:** tests and scripts run inside the container:

```bash
docker run --rm --pull always -v "$PWD":/work -w /work \
  veupathdb/dnaseqanalysis:latest bash -c 'bash testing/t/singleSampleMerge.t.sh'
```

---

## File Structure

| File | Responsibility |
|---|---|
| Create `bin/mergeVcfs.sh` | Produce the `merged.vcf.gz` contract from 1..N per-strain VCFs; error on 0 |
| Create `bin/mergeCoverageBeds.sh` | Produce the `coverage.tsv` contract from 1..N per-strain coverage BEDs; error on 0 |
| Modify `modules/mergeExperiments.nf` | `mergeVcfs` and `mergeCoverageBeds` process bodies become one-line script calls |
| Modify `workflows/mergeExperiments.nf` | Delete the `single`/`multiple` branch |
| Create `testing/t/singleSampleMerge.t.sh` | Characterization tests for both scripts at n=0, n=1, n=2 |
| Create `testing/bin/compareMergeOutputs.py` | Compare two output dirs modulo strain-id permutation (an instrument, not a test — lives outside `testing/t/` so pytest does not collect it) |

---

## Task 1: Extract `mergeVcfs` shell body into `bin/mergeVcfs.sh`

Pure extraction. No behavior change yet — n=1 still fails after this task. Doing
the move first means the next task's test has something to call.

**Files:**
- Create: `bin/mergeVcfs.sh`
- Modify: `modules/mergeExperiments.nf:26-51`

- [ ] **Step 1: Create the script with today's logic verbatim**

```bash
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
```

- [ ] **Step 2: Make it executable**

```bash
chmod +x bin/mergeVcfs.sh
```

- [ ] **Step 3: Replace the process body**

In `modules/mergeExperiments.nf`, replace the whole `script:` block of
`mergeVcfs` (currently lines 35-45) with:

```groovy
  script:
    """
    set -euo pipefail
    mergeVcfs.sh merged.vcf.gz *.vcf.gz
    """
```

Leave `container`, `input`, `output` and `stub` exactly as they are. The glob is
expanded by the calling shell before `mergeVcfs.sh` runs, so the `.norm.vcf.gz`
files the script creates are never picked up as inputs.

- [ ] **Step 4: Verify the extraction preserved behavior on 2 samples**

```bash
docker run --rm -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest bash -c '
set -e
cd "$(mktemp -d)"
hdr() { printf "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n##INFO=<ID=AO,Number=A,Type=Integer,Description=\"AO\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" "$1"; }
{ hdr SA; printf "chr1\t100\t.\tA\tG\t.\t.\tAO=5\tGT\t1/1\n"; } > 1.vcf
{ hdr SB; printf "chr1\t200\t.\tC\tT\t.\t.\tAO=7\tGT\t0/1\n"; } > 2.vcf
bgzip 1.vcf; bgzip 2.vcf
/work/bin/mergeVcfs.sh merged.vcf.gz *.vcf.gz
bcftools view -h merged.vcf.gz | tail -1
bcftools view -H merged.vcf.gz
'
```

Expected: two records, sample columns `SA` and `SB` present, `INFO` column `.`
on both rows.

- [ ] **Step 5: Commit**

```bash
git add bin/mergeVcfs.sh modules/mergeExperiments.nf
git commit -m "Extract mergeVcfs shell body into bin/mergeVcfs.sh"
```

---

## Task 2: `bin/mergeVcfs.sh` handles a single input

**Files:**
- Create: `testing/t/singleSampleMerge.t.sh`
- Modify: `bin/mergeVcfs.sh`

- [ ] **Step 1: Write the failing test file**

Create `testing/t/singleSampleMerge.t.sh`. This is the whole file; later tasks
append to it.

```bash
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
```

- [ ] **Step 2: Run it and watch it fail**

```bash
docker run --rm -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest \
  bash testing/t/singleSampleMerge.t.sh
```

Expected: `FAIL: mergeVcfs.sh exited non-zero on a single input` (the underlying
cause is `bcftools merge`'s usage dump — it needs two or more files).

- [ ] **Step 3: Make the combine step arity-aware**

In `bin/mergeVcfs.sh`, replace the final line

```bash
bcftools merge --merge both -O z -o "$out" "${norm[@]}"
```

with

```bash
if (( ${#norm[@]} == 1 )); then
  # Single strain: nothing to merge, and bcftools merge rejects one input. The
  # normalized file already IS the contract -- verified that bcftools merge does
  # not re-collapse rows that norm -m -any split within a single file, so this
  # is shape-identical to the n>=2 output minus the other samples' columns.
  mv "${norm[0]}" "$out"
else
  bcftools merge --merge both -O z -o "$out" "${norm[@]}"
fi
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
docker run --rm -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest \
  bash testing/t/singleSampleMerge.t.sh
```

Expected: `PASS: mergeVcfs.sh n=1`

- [ ] **Step 5: Commit**

```bash
git add bin/mergeVcfs.sh testing/t/singleSampleMerge.t.sh
git commit -m "mergeVcfs.sh: emit the normalized VCF directly when there is one strain"
```

---

## Task 3: `bin/mergeVcfs.sh` still merges 2+ inputs, and rejects 0

Guards against a fix that quietly breaks the multi-sample path, and adds the
fail-fast on an empty input glob.

**Files:**
- Modify: `testing/t/singleSampleMerge.t.sh`
- Modify: `bin/mergeVcfs.sh`

- [ ] **Step 1: Append the failing tests**

Append to `testing/t/singleSampleMerge.t.sh`:

```bash
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
```

- [ ] **Step 2: Run and confirm which case fails**

```bash
docker run --rm -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest \
  bash testing/t/singleSampleMerge.t.sh
```

Expected: n=1 and n=2 pass; the n=0 case fails. With `norm` empty,
`"${norm[@]}"` under `set -u` makes the script die on an unbound-variable error
rather than a message naming the problem, so expect
`FAIL: mergeVcfs.sh n=0 message did not name the missing input: ...`.

- [ ] **Step 3: Add the n=0 guard**

In `bin/mergeVcfs.sh`, immediately after the `shift` that consumes `$out`, insert:

```bash
if (( $# == 0 )); then
  echo "ERROR: mergeVcfs.sh: no input VCFs -- mergeExperiments needs at least one per-strain VCF (check params.vcfFiles)" >&2
  exit 1
fi
```

- [ ] **Step 4: Run the test to verify all three cases pass**

```bash
docker run --rm -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest \
  bash testing/t/singleSampleMerge.t.sh
```

Expected: `PASS: mergeVcfs.sh n=1`, `PASS: mergeVcfs.sh n=2`,
`PASS: mergeVcfs.sh n=0`.

- [ ] **Step 5: Commit**

```bash
git add bin/mergeVcfs.sh testing/t/singleSampleMerge.t.sh
git commit -m "mergeVcfs.sh: fail loudly when no input VCFs were staged"
```

---

## Task 4: Extract `mergeCoverageBeds` shell body into `bin/mergeCoverageBeds.sh`

Pure extraction, same shape as Task 1.

**Files:**
- Create: `bin/mergeCoverageBeds.sh`
- Modify: `modules/mergeExperiments.nf:55-81`

- [ ] **Step 1: Create the script**

```bash
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
```

Note on the header: the old process body wrote it via `IFS='\t'` (a literal
backslash and `t`, not a tab) and relied on `echo -e` to turn the resulting
`A\tB` back into tabs. `IFS=$'\t'` plus `printf` produces the same bytes
directly. Output is unchanged; the accident is just no longer load-bearing.

- [ ] **Step 2: Make it executable**

```bash
chmod +x bin/mergeCoverageBeds.sh
```

- [ ] **Step 3: Replace the process body**

In `modules/mergeExperiments.nf`, replace the whole `script:` block of
`mergeCoverageBeds` (currently lines 64-75) with:

```groovy
  script:
    """
    set -euo pipefail
    mergeCoverageBeds.sh coverage.tsv $coverageBeds
    """
```

`$coverageBeds` is the Nextflow input variable and interpolates to the staged
filenames — leave it unescaped. `container`, `input`, `output` and `stub` are
unchanged.

- [ ] **Step 4: Verify the extraction preserved behavior on 2 samples**

```bash
docker run --rm -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest bash -c '
set -e
cd "$(mktemp -d)"
printf "chr1\t0\t100\t5\n"   | bgzip > StrainA_coverage.bed.gz
printf "chr1\t50\t150\t9\n"  | bgzip > StrainB_coverage.bed.gz
/work/bin/mergeCoverageBeds.sh coverage.tsv StrainA_coverage.bed.gz StrainB_coverage.bed.gz
cat -A coverage.tsv | head -1
cat coverage.tsv
'
```

Expected: header line is `chrom^Istart^Iend^IStrainA^IStrainB$` under `cat -A`
(real tabs, no literal `\t`), followed by three union intervals.

- [ ] **Step 5: Commit**

```bash
git add bin/mergeCoverageBeds.sh modules/mergeExperiments.nf
git commit -m "Extract mergeCoverageBeds shell body into bin/mergeCoverageBeds.sh"
```

---

## Task 5: `bin/mergeCoverageBeds.sh` handles a single input, and rejects 0

**Files:**
- Modify: `testing/t/singleSampleMerge.t.sh`
- Modify: `bin/mergeCoverageBeds.sh`

- [ ] **Step 1: Append the failing tests**

Append to `testing/t/singleSampleMerge.t.sh`:

```bash
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
```

- [ ] **Step 2: Run and watch the n=1 and n=0 cases fail**

```bash
docker run --rm -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest \
  bash testing/t/singleSampleMerge.t.sh
```

Expected: mergeVcfs cases pass, `PASS: mergeCoverageBeds.sh n=2` passes, and
`FAIL: mergeCoverageBeds.sh exited non-zero on a single input` (from
`bedtools unionbedg`'s "Only a single BedGraph file was specified").

- [ ] **Step 3: Add the n=0 guard and the n=1 branch**

In `bin/mergeCoverageBeds.sh`, immediately after the `shift`, insert:

```bash
if (( $# == 0 )); then
  echo "ERROR: mergeCoverageBeds.sh: no input coverage BEDs -- mergeExperiments needs at least one per-strain coverage BED (check params.coverageFiles)" >&2
  exit 1
fi
```

Then replace the final line

```bash
bedtools unionbedg -names "${names[@]}" -filler 0 -i "$@" >> "$out"
```

with

```bash
if (( $# == 1 )); then
  # Single strain: nothing to union, and unionbedg rejects one input. Its output
  # for a single 4-column bed is byte-identical to that file's own contents.
  zcat "$1" >> "$out"
else
  bedtools unionbedg -names "${names[@]}" -filler 0 -i "$@" >> "$out"
fi
```

- [ ] **Step 4: Run the whole suite**

```bash
docker run --rm -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest \
  bash testing/t/singleSampleMerge.t.sh
```

Expected, in order: `PASS: mergeVcfs.sh n=1`, `PASS: mergeVcfs.sh n=2`,
`PASS: mergeVcfs.sh n=0`, `PASS: mergeCoverageBeds.sh n=1`,
`PASS: mergeCoverageBeds.sh n=2`, `PASS: mergeCoverageBeds.sh n=0`.

- [ ] **Step 5: Commit**

```bash
git add bin/mergeCoverageBeds.sh testing/t/singleSampleMerge.t.sh
git commit -m "mergeCoverageBeds.sh: handle a single strain, fail loudly on none"
```

---

## Task 6: Delete the workflow-level arity branch

This is the task that removes the wrong-VCF bug. Until now the single-sample
bypass was still in place and would still have shadowed the fixed script.

**Files:**
- Modify: `workflows/mergeExperiments.nf:26-33`

- [ ] **Step 1: Replace the branch with a direct call**

In `workflows/mergeExperiments.nf`, replace these lines

```groovy
    allVcfs       = vcfs_qch.collect()
    allVcfsBranch = allVcfs.branch { single: it.size() == 1; multiple: true }
    mergedVcf     = allVcfsBranch.single.map { it[0] }
                      .mix(mergeVcfs(allVcfsBranch.multiple))
```

with

```groovy
    // Arity is mergeVcfs' problem, not the DAG's: bin/mergeVcfs.sh normalizes
    // every input and merges only when there are 2+. Branching here previously
    // let the n=1 case skip normalization entirely.
    mergedVcf = mergeVcfs(vcfs_qch.collect())
```

Leave `allFastas`, `combinedIndels`, and the `mergeCoverageBeds` call as they
are — `mergeCoverageBeds(coverages_qch.collect())` already has the right shape.

- [ ] **Step 2: Verify both entry points still parse**

```bash
nextflow run main.nf -entry mergeExperiments -profile mergeExperiments -stub-run \
  -c /home/jbrestel/dnaseq_test/single_sample/nextflow.config 2>&1 | tail -20
```

Expected: the DAG resolves and the stub run completes. A `Missing process or
function` / `Unknown method` error means the edit is wrong. (Stub bodies write
empty files, so output content is meaningless here — this checks wiring only.)

- [ ] **Step 3: Commit**

```bash
git add workflows/mergeExperiments.nf
git commit -m "Drop the single/multiple VCF branch; mergeVcfs owns its arity"
```

---

## Task 7: Write the strain-id-permutation-aware comparator

Needed because strain ids are assigned in channel-staging order, which is not
stable across runs — a plain `diff` against a baseline reports differences that
are pure renumbering. See the "Strain ids are run-order dependent" section of
the spec.

Two distinct numberings exist and must each be read from their own file:
`sample.dat` (`strain_id`→`sample_name`, reference last) and
`hsss_readFreq*/strainIdToName.dat` (reference first).

**Files:**
- Create: `testing/bin/compareMergeOutputs.py`

- [ ] **Step 1: Create the comparator**

```python
#!/usr/bin/env python3
"""Compare two mergeExperiments output directories for equivalence modulo
strain-id permutation.

Strain ids are assigned in whatever order the collected per-strain files were
staged, so two runs over identical inputs can produce identical data under
different numbering. This canonicalizes every id-bearing field to sample NAMES
before comparing, and sorts rows (row order within a file is not guaranteed
either).

Usage: compareMergeOutputs.py DIR_A DIR_B
Exit 0 if equivalent, 1 if not, 2 on a usage/IO problem.
"""

import subprocess
import sys
from pathlib import Path

# .dat files with no id-bearing column: compared as-is (sorted).
PLAIN_DATS = ["variationFeature.dat", "snpeff.dat"]

# .dat files with a column holding a "{1,3}"-style strain-id set.
ID_SET_DATS = {
    "allele.dat": "strain_ids",
    "transcript_product.dat": "downstream_of_frameshift_strain_ids",
}

problems = []


def problem(msg):
    problems.append(msg)


def read_id_map(path):
    """Read a two-column id<TAB>name file into {id: name}, skipping a header
    row if present. Used for both sample.dat and strainIdToName.dat."""
    mapping = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            problem(f"{path}: cannot parse line: {line!r}")
            continue
        sid, name = parts[0], parts[1]
        if sid == "strain_id":  # sample.dat header
            continue
        mapping[sid] = name
    return mapping


def canon_id_set(field, id_map, where):
    """'{1,3}' -> '{Friedlin,LV39}' (name-sorted). Empty stays empty."""
    raw = field.strip()
    if not raw or raw == ".":
        return raw
    if not (raw.startswith("{") and raw.endswith("}")):
        problem(f"{where}: expected a {{...}} id set, got {field!r}")
        return raw
    ids = [i for i in raw[1:-1].split(",") if i]
    names = []
    for i in ids:
        if i not in id_map:
            problem(f"{where}: strain id {i} not in sample.dat")
            names.append(f"<unknown:{i}>")
        else:
            names.append(id_map[i])
    return "{" + ",".join(sorted(names)) + "}"


def load_dat(path, id_col=None, id_map=None):
    """Return (header, sorted canonicalized data rows)."""
    lines = path.read_text().splitlines()
    if not lines:
        return "", []
    header, rows = lines[0], lines[1:]
    if id_col is None:
        return header, sorted(rows)
    cols = header.lstrip("#").split("\t")
    if id_col not in cols:
        problem(f"{path}: no {id_col} column in header")
        return header, sorted(rows)
    idx = cols.index(id_col)
    out = []
    for row in rows:
        fields = row.split("\t")
        if idx < len(fields):
            fields[idx] = canon_id_set(fields[idx], id_map, f"{path.name}")
        out.append("\t".join(fields))
    return header, sorted(out)


def compare_rows(name, a, b):
    ah, arows = a
    bh, brows = b
    if ah != bh:
        problem(f"{name}: header differs\n  A: {ah}\n  B: {bh}")
    if arows == brows:
        return
    only_a = [r for r in arows if r not in set(brows)]
    only_b = [r for r in brows if r not in set(arows)]
    problem(
        f"{name}: {len(only_a)} row(s) only in A, {len(only_b)} row(s) only in B\n"
        + "".join(f"  A only: {r}\n" for r in only_a[:5])
        + "".join(f"  B only: {r}\n" for r in only_b[:5])
    )


def compare_sample_dat(a_dir, b_dir):
    """Sample identity is the set of names; the id assignment is what we are
    deliberately ignoring."""
    a = set(read_id_map(a_dir / "sample.dat").values())
    b = set(read_id_map(b_dir / "sample.dat").values())
    if a != b:
        problem(f"sample.dat: sample sets differ\n  only A: {sorted(a - b)}\n  only B: {sorted(b - a)}")


def compare_hsss(a_dir, b_dir):
    """hsss dirs hold per-strain files NAMED by a strain index that uses its own
    numbering (reference first). Map filenames to names via each dir's own
    strainIdToName.dat, then compare contents pairwise by name."""
    for freq in ("20", "40", "60", "80"):
        d = f"hsss_readFreq{freq}"
        a, b = a_dir / d, b_dir / d
        if not a.is_dir() or not b.is_dir():
            problem(f"{d}: missing in {'A' if not a.is_dir() else 'B'}")
            continue
        for meta in ("contigIdToSourceId.dat", "referenceGenome.dat"):
            if (a / meta).read_bytes() != (b / meta).read_bytes():
                problem(f"{d}/{meta}: differs")
        a_map = read_id_map(a / "strainIdToName.dat")
        b_map = read_id_map(b / "strainIdToName.dat")
        if set(a_map.values()) != set(b_map.values()):
            problem(f"{d}/strainIdToName.dat: strain sets differ")
            continue
        b_by_name = {name: sid for sid, name in b_map.items()}
        for sid, name in a_map.items():
            af, bf = a / sid, b / b_by_name[name]
            if not af.exists() or not bf.exists():
                problem(f"{d}: no per-strain file for {name} (A:{af.name} B:{bf.name})")
                continue
            if af.read_bytes() != bf.read_bytes():
                problem(f"{d}: per-strain data for {name} differs ({af.name} vs {bf.name})")


def vcf_body(path, sample_order):
    """Record body only, samples forced into a canonical order. The ## header
    carries bcftools/snpEff version and command-line stamps, so the compressed
    files differ even when the data does not."""
    return subprocess.run(
        ["bcftools", "view", "-H", "-s", ",".join(sample_order), str(path)],
        capture_output=True, text=True, check=True,
    ).stdout


def compare_ann_vcf(a_dir, b_dir):
    name = "merged.ann.vcf.gz"
    a, b = a_dir / name, b_dir / name
    if not a.exists() or not b.exists():
        problem(f"{name}: missing in {'A' if not a.exists() else 'B'}")
        return

    def samples(p):
        out = subprocess.run(["bcftools", "query", "-l", str(p)],
                             capture_output=True, text=True, check=True)
        return sorted(out.stdout.split())

    sa, sb = samples(a), samples(b)
    if sa != sb:
        problem(f"{name}: sample sets differ\n  A: {sa}\n  B: {sb}")
        return
    if sorted(vcf_body(a, sa).splitlines()) != sorted(vcf_body(b, sb).splitlines()):
        problem(f"{name}: record body differs")


def main():
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    a_dir, b_dir = Path(sys.argv[1]), Path(sys.argv[2])
    for d in (a_dir, b_dir):
        if not d.is_dir():
            print(f"not a directory: {d}", file=sys.stderr)
            return 2

    a_map = read_id_map(a_dir / "sample.dat")
    b_map = read_id_map(b_dir / "sample.dat")

    compare_sample_dat(a_dir, b_dir)

    for name in PLAIN_DATS:
        compare_rows(name, load_dat(a_dir / name), load_dat(b_dir / name))

    for name, col in ID_SET_DATS.items():
        compare_rows(name,
                     load_dat(a_dir / name, col, a_map),
                     load_dat(b_dir / name, col, b_map))

    compare_hsss(a_dir, b_dir)
    compare_ann_vcf(a_dir, b_dir)

    if problems:
        print(f"DIFFER: {len(problems)} problem(s)\n")
        for p in problems:
            print(f"- {p}")
        return 1
    print("EQUIVALENT (modulo strain-id permutation)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2: Make it executable**

```bash
mkdir -p testing/bin && chmod +x testing/bin/compareMergeOutputs.py
```

- [ ] **Step 3: Check it calls a directory equivalent to itself**

```bash
docker run --rm -v "$PWD":/work -v /home/jbrestel/dnaseq_test:/data -w /work \
  veupathdb/dnaseqanalysis:latest \
  python3 testing/bin/compareMergeOutputs.py /data/merge/output.orig /data/merge/output.orig
```

Expected: `EQUIVALENT (modulo strain-id permutation)`, exit 0. If this reports
differences, the comparator is wrong — a directory is always equivalent to
itself.

- [ ] **Step 4: Check it detects a real difference**

```bash
docker run --rm -v "$PWD":/work -v /home/jbrestel/dnaseq_test:/data -w /work \
  veupathdb/dnaseqanalysis:latest bash -c '
set -e
t="$(mktemp -d)"; cp -r /data/merge/output.orig "$t/edited"
# perturb one allele frequency
awk -F"\t" -v OFS="\t" "NR==2 { \$5 = \"0.1234\" } 1" \
  /data/merge/output.orig/allele.dat > "$t/edited/allele.dat"
python3 testing/bin/compareMergeOutputs.py /data/merge/output.orig "$t/edited"
'
```

Expected: `DIFFER: 1 problem(s)` naming `allele.dat`, exit 1. A comparator that
passes here is useless.

- [ ] **Step 5: Check it is blind to a consistent id permutation**

This is the whole point of the comparator, so it gets its own fixture. The
script below rewrites a copy of the baseline with strain ids 1 and 3 swapped
consistently across all four id-bearing artifacts, in both numbering schemes.

```bash
docker run --rm -v "$PWD":/work -v /home/jbrestel/dnaseq_test:/data -w /work \
  veupathdb/dnaseqanalysis:latest python3 - <<'PY'
import re, shutil, subprocess, sys, tempfile
from pathlib import Path

src = Path("/data/merge/output.orig")
dst = Path(tempfile.mkdtemp()) / "permuted"
shutil.copytree(src, dst)

def swap_ids(text, a="1", b="3"):
    """Swap two ids inside every {..} set in the text."""
    def repl(m):
        ids = [b if i == a else a if i == b else i for i in m.group(1).split(",") if i]
        return "{" + ",".join(ids) + "}"
    return re.sub(r"\{([0-9,]*)\}", repl, text)

# sample.dat: swap which name each id points at
lines = (dst / "sample.dat").read_text().splitlines()
out = [lines[0]]
for line in lines[1:]:
    sid, name = line.split("\t")
    out.append("\t".join(("3" if sid == "1" else "1" if sid == "3" else sid, name)))
(dst / "sample.dat").write_text("\n".join(out) + "\n")

# allele.dat + transcript_product.dat: swap ids inside the {..} sets
for f in ("allele.dat", "transcript_product.dat"):
    p = dst / f
    p.write_text(swap_ids(p.read_text()))

# hsss dirs: swap the per-strain FILE names and the strainIdToName rows.
# This numbering is independent of sample.dat's, so swapping 1<->3 here is a
# second, unrelated permutation -- which the comparator must also absorb.
for freq in ("20", "40", "60", "80"):
    d = dst / f"hsss_readFreq{freq}"
    (d / "1").rename(d / "tmp"); (d / "3").rename(d / "1"); (d / "tmp").rename(d / "3")
    rows = []
    for line in (d / "strainIdToName.dat").read_text().splitlines():
        sid, name = line.split("\t")
        rows.append("\t".join(("3" if sid == "1" else "1" if sid == "3" else sid, name)))
    (d / "strainIdToName.dat").write_text("\n".join(rows) + "\n")

rc = subprocess.run([sys.executable, "testing/bin/compareMergeOutputs.py",
                     str(src), str(dst)]).returncode
sys.exit(rc)
PY
```

Expected: `EQUIVALENT (modulo strain-id permutation)`, exit 0. If it reports
differences, the comparator is still id-order sensitive and the Task 8
regression check would give false alarms.

- [ ] **Step 6: Commit**

```bash
git add testing/bin/compareMergeOutputs.py
git commit -m "Add strain-id-permutation-aware mergeExperiments output comparator"
```

---

## Task 8: Multi-sample regression — prove nothing changed at n>=2

The n>=2 path is untouched by design, so this must show zero change. Any
difference means one of Tasks 1-6 altered multi-sample behavior and is wrong.

**Files:** none modified. This is a verification task.

- [ ] **Step 1: Rerun the 4-sample dataset**

Do **not** pass `-resume` — cached task results would defeat the point.

```bash
cd /home/jbrestel/dnaseq_test/merge
nextflow -C $PWD/nextflow.config run \
  /home/jbrestel/workspaces/dataLoad/project_home/dnaseq-nextflow/main.nf \
  -entry mergeExperiments -profile mergeExperiments
```

Expected: the run completes; `nextflow log | tail -n 1 | cut -f 3,4` reports
status `OK`. Output lands in `/home/jbrestel/dnaseq_test/merge/output/`.

- [ ] **Step 2: Compare against the baseline**

```bash
cd /home/jbrestel/workspaces/dataLoad/project_home/dnaseq-nextflow
docker run --rm -v "$PWD":/work -v /home/jbrestel/dnaseq_test:/data -w /work \
  veupathdb/dnaseqanalysis:latest \
  python3 testing/bin/compareMergeOutputs.py /data/merge/output.orig /data/merge/output
```

Expected: `EQUIVALENT (modulo strain-id permutation)`, exit 0.

If it reports differences: stop and investigate before continuing. The likely
culprits are the header-construction change in Task 4 (compare
`work/*/coverage.tsv` between the old and new runs with `cat -A`) or the
`mergeVcfs.sh` argument order changing the sample column order — the latter is
a permutation the comparator absorbs, so it should not surface here at all.

- [ ] **Step 3: Record the result**

No commit — nothing changed. Note in the task handoff whether the comparator
reported `EQUIVALENT`, and paste its actual output. Do not claim this passed
without the command's output in hand.

---

## Task 9: Single-sample end-to-end run

The reason this whole plan exists. Until now nothing has proven the pipeline runs
end to end with one sample.

**Files:** none modified. This is a verification task.

- [ ] **Step 1: Run the single-sample dataset**

```bash
cd /home/jbrestel/dnaseq_test/single_sample
nextflow -C $PWD/nextflow.config run \
  /home/jbrestel/workspaces/dataLoad/project_home/dnaseq-nextflow/main.nf \
  -entry mergeExperiments -profile mergeExperiments
```

Expected: completes through `parseSnpEffAnnotations`. Verify with:

```bash
nextflow log | tail -n 1 | cut -f 3,4
```

Expected: status `OK`.

If a process failed, get its work directory and read `.command.err`:

```bash
runName=$(nextflow log | tail -n 1 | cut -f 3)
nextflow log $runName -f native_id,exit,process,workdir | grep -v ' 0 '
```

- [ ] **Step 2: Confirm the merged VCF got the contract, not the raw bypass**

This is the specific bug the plan fixes, so check it directly on the real data.

```bash
cd /home/jbrestel/dnaseq_test/single_sample
runName=$(nextflow log | tail -n 1 | cut -f 3)
wd=$(nextflow log $runName -f process,workdir | grep mergeVcfs | cut -f2)
docker run --rm -v "$wd":/wd veupathdb/dnaseqanalysis:latest bash -c '
bcftools view -H /wd/merged.vcf.gz | cut -f8 | sort -u | head -3
echo "--- FORMAT keys:"
bcftools view -H /wd/merged.vcf.gz | cut -f9 | sort -u
echo "--- multiallelic ALT rows (want 0):"
bcftools view -H /wd/merged.vcf.gz | cut -f5 | grep -c "," || true
'
```

Expected: INFO column is `.` on every row, FORMAT keys contain neither `GL` nor
`DPR`, and zero multiallelic ALT rows. Anything else means the normalization
still is not reaching the single-sample path.

- [ ] **Step 3: Confirm the published artifacts are sane for one strain**

```bash
cd /home/jbrestel/dnaseq_test/single_sample/output
cat sample.dat
echo "--- rows per file:"
wc -l allele.dat variationFeature.dat transcript_product.dat snpeff.dat
echo "--- hsss:"
ls hsss_readFreq20 hsss_readFreq40 hsss_readFreq60 hsss_readFreq80
echo "--- allele.dat strain_ids seen:"
cut -f8 allele.dat | sort -u
```

Expected:
- `sample.dat` has exactly two rows: `Seidman751` and the reference
  `lmajFriedlin`. (`reference_strain` is external —
  `bin/processSequenceVariations.jl` appends it to `all_strains` rather than
  selecting it from them — so one strain plus the reference is well-formed.)
- All four `.dat` files are non-empty (more than just a header line).
- Each `hsss_readFreq*` directory contains `contigIdToSourceId.dat`,
  `referenceGenome.dat`, `strainIdToName.dat`, and per-strain files.
- `strain_ids` sets reference only the two ids from `sample.dat`.

- [ ] **Step 4: Sanity-read Seidman751 against the 4-sample run**

Not a diff — the two runs legitimately differ. Frequency and count columns use a
distinct-strain denominator (1 strain vs 4), so `allele_frequency`,
`distinct_strain_count` and the `*_strain_count` columns *should* differ. What
should agree is the per-strain evidence: at a position where Seidman751 carries a
non-reference allele in the 4-sample run, the same allele should appear in the
single-sample run with a similar `avg_coverage` / `avg_percent`.

```bash
# pick a position where Seidman751 (id 3 in the baseline) is the only alt carrier
grep -P '\t\{3\}\t' /home/jbrestel/dnaseq_test/merge/output.orig/allele.dat | head -5
# then look up those locations in the single-sample run
```

Expected: the allele letter, `avg_coverage` and `avg_percent` for those locations
are close to the 4-sample values (coverage is per-strain, so it should be
essentially identical). A wildly different coverage figure means `coverage.tsv`
is being misread at n=1 — investigate the sample-name column header in the
single-sample `coverage.tsv`.

- [ ] **Step 5: Report findings**

If `bin/processSequenceVariations.jl` itself misbehaves at n=1 (degenerate
major/minor allele calls, frequency edge cases), **report it and stop** — per the
spec, fixing the Julia code is out of scope and gets its own spec. Do not patch
it here.

---

## Task 10: Document the new suite and finish the branch

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Add the comparator to the testing docs**

In `CLAUDE.md`, in the `## Testing` section, after the paragraph describing the
`--run-dir` e2e tests, add:

```markdown
`testing/bin/compareMergeOutputs.py DIR_A DIR_B` compares two `mergeExperiments`
output directories for equivalence modulo strain-id permutation. Strain ids are
assigned in channel-staging order, so two runs over identical inputs can produce
identical data under different numbering — a plain `diff` reports differences
that are pure renumbering. Use this instead when checking that a change left
multi-sample output untouched.
```

The `for t in testing/t/*.t.sh` loop already documented in `CLAUDE.md` picks up
`singleSampleMerge.t.sh` with no change needed.

- [ ] **Step 2: Run the full test suite one last time**

```bash
docker run --rm --pull always -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest bash -c '
  for t in testing/t/*.jl; do julia "$t"; done
  python3 -m pytest testing/t/
  for t in testing/t/*.t.sh; do echo "== $t"; bash "$t"; done
'
```

Expected: all Julia suites pass, pytest passes (the 96 `test_mergeExperiments_e2e.py`
tests skip without `--run-dir`), and all three `*.t.sh` suites print `PASS`
lines. Paste the actual output into the handoff — do not assert this passed
without it.

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md
git commit -m "Document the mergeExperiments output comparator"
```

- [ ] **Step 4: Finish the branch**

Use the superpowers:finishing-a-development-branch skill to decide between merge,
PR, or further work. Do not merge to `main` without asking.

---

## Verification summary

What must be true before calling this done, with the evidence for each:

| Claim | Evidence |
|---|---|
| n=1 satisfies the same VCF contract as n>=2 | `testing/t/singleSampleMerge.t.sh` passes (Tasks 2-3) |
| n=1 satisfies the same coverage contract | same suite (Task 5) |
| n=0 fails loudly instead of confusingly | same suite (Tasks 3, 5) |
| The wrong-VCF bypass is gone | Task 6 edit + Task 9 Step 2 on real data |
| n>=2 behavior is unchanged | Task 8 comparator output says `EQUIVALENT` |
| A single-sample run completes end to end | Task 9 Step 1 reports status `OK` |
| Nothing else in the suite regressed | Task 10 Step 2 output |
