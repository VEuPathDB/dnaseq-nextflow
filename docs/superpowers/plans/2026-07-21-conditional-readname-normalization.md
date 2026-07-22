# Conditional Read-Name Normalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Only run FASTQ read-name normalization (and its lockstep guard) on paired samples whose read names carry a mate suffix `bwa-mem2` cannot handle; pass already-conforming samples to bwa untouched.

**Architecture:** Extend `bin/normalizeReadNames.sh` with a `needs-normalizing` predicate mode (single source of truth for the suffix set). In `modules/alignment.nf`, the paired branch full-scans both mates with that predicate; only when it fires does it run the lockstep length assertion and feed mates through the `strip` filter via process substitution — otherwise it hands raw files straight to bwa.

**Tech Stack:** Bash / awk, Nextflow DSL2, `veupathdb/dnaseqanalysis:latest` container for tests.

**Working directory:** `.worktrees/fix-mate-readname-mismatch` (branch `fix-mate-readname-mismatch`, extends PR #25).

**Spec:** `docs/superpowers/specs/2026-07-21-conditional-readname-normalization-design.md`

---

## Task 1: Add `needs-normalizing` detection mode to `normalizeReadNames.sh`

**Files:**
- Modify: `bin/normalizeReadNames.sh`
- Create: `testing/t/normalizeReadNames.t.sh`

- [ ] **Step 1: Write the failing test**

Create `testing/t/normalizeReadNames.t.sh`:

```bash
#!/usr/bin/env bash
# Unit tests for bin/normalizeReadNames.sh: the `strip` transform and the
# `needs-normalizing` predicate. No -e: `needs-normalizing` uses exit codes as
# data, and command-substitution pipelines below are inspected, not asserted.
set -uo pipefail

SCRIPT="$(cd "$(dirname "$0")/../.." && pwd)/bin/normalizeReadNames.sh"
fail() { echo "FAIL: $1"; exit 1; }

# a FASTQ record is 4 lines: @id, seq, +, qual. build one with a given id.
rec() { printf '%s\nACGT\n+\nFFFF\n' "$1"; }
# variant whose QUALITY line ends in ":1" (a suffix-looking string on a non-ID line)
recq() { printf '%s\nACGT\n+\nFF:1\n' "$1"; }

# --- strip mode --------------------------------------------------------------
out="$(rec '@read/1' | "$SCRIPT" strip)";        [ "${out%%$'\n'*}" = "@read" ]          || fail "strip /1: '${out%%$'\n'*}'"
out="$(rec '@SRR868698.4.1' | "$SCRIPT")";       [ "${out%%$'\n'*}" = "@SRR868698.4" ]   || fail "strip .1 (default mode): '${out%%$'\n'*}'"
out="$(rec '@90206:7:1:1:518:a' | "$SCRIPT" strip)"; [ "${out%%$'\n'*}" = "@90206:7:1:1:518" ] || fail "strip :a: '${out%%$'\n'*}'"
out="$(rec '@keep/1extra' | "$SCRIPT" strip)";   [ "${out%%$'\n'*}" = "@keep/1extra" ]   || fail "strip keep/1extra: '${out%%$'\n'*}'"
out="$(rec '@plain.name' | "$SCRIPT" strip)";    [ "${out%%$'\n'*}" = "@plain.name" ]    || fail "strip plain.name: '${out%%$'\n'*}'"

# only every 4th line is touched: id stripped, qual ":1" left intact
out="$(recq '@read.1' | "$SCRIPT" strip)"
[ "${out%%$'\n'*}" = "@read" ]  || fail "strip recq id: '${out%%$'\n'*}'"
[ "${out##*$'\n'}" = "FF:1" ]   || fail "strip touched a non-ID line: '${out##*$'\n'}'"

# --- needs-normalizing mode (exit 0 = needs, exit 1 = conforming) ------------
if   { rec '@read/1'; rec '@read/2'; } | "$SCRIPT" needs-normalizing; then fail "/1,/2 should be conforming (exit 1)"; fi
if   { rec '@read';   rec '@read';   } | "$SCRIPT" needs-normalizing; then fail "identical names should be conforming (exit 1)"; fi
if ! { rec '@x.1';    rec '@x.2';    } | "$SCRIPT" needs-normalizing; then fail ".1,.2 should need normalizing (exit 0)"; fi
if ! { rec '@x:a';    rec '@x:b';    } | "$SCRIPT" needs-normalizing; then fail ":a,:b should need normalizing (exit 0)"; fi

# mixed file: conforming identical-name reads (end in a coordinate, not a mate
# suffix) followed by one :a/:b pair -> whole sample needs normalizing (exit 0)
if ! { rec '@090126:3:1:1:1092'; rec '@090126:3:1:1:1092'; \
       rec '@90206:7:1:1:518:a'; rec '@90206:7:1:1:518:b'; } \
     | "$SCRIPT" needs-normalizing; then fail "mixed file should need normalizing (exit 0)"; fi

# --- usage -------------------------------------------------------------------
rc=0; "$SCRIPT" bogus </dev/null 2>/dev/null || rc=$?
[ "$rc" -eq 2 ] || fail "unknown mode should exit 2, got $rc"

echo "PASS"
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
docker run --rm --pull always -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest \
  bash testing/t/normalizeReadNames.t.sh
```
Expected: FAIL — the current script ignores its argument and always runs the strip
awk, so `"$SCRIPT" needs-normalizing` strips instead of predicating (exits 0 for a
conforming `/1,/2` file), and `"$SCRIPT" bogus` exits 0 instead of 2. First
failing assertion: `/1,/2 should be conforming (exit 1)`.

- [ ] **Step 3: Rewrite `bin/normalizeReadNames.sh` with mode dispatch**

Replace the entire file with:

```bash
#!/bin/bash
#
# normalizeReadNames.sh - strip trailing mate suffixes from FASTQ read-ID lines,
# or detect whether a FASTQ carries any such suffix.
#
# bwa-mem2 pairs mate1/mate2 by requiring IDENTICAL read names; it only auto-strips
# a trailing /1 or /2. Reads carrying other mate suffixes therefore abort alignment
# with "[mem_sam_pe*] paired reads have different names". Two conventions show up in
# our source data and both trip it:
#   - SRA fastq-dump --readids  ->  .1 / .2   (e.g. SRR868698.4.1 / SRR868698.4.2)
#   - legacy Illumina/Casava    ->  :a / :b   (e.g. 90206:7:1:1:518:a / :b)
#
# Modes (first argument, default "strip"):
#   strip              read FASTQ on stdin, strip a trailing /1 /2, .1 .2, or :a :b
#                      from every read-ID line (every 4th line), write to stdout.
#                      Pure stream filter (no temp files) so it can run inside a
#                      process substitution feeding bwa-mem2 at zero extra I/O cost.
#   needs-normalizing  read FASTQ on stdin, exit 0 if any read-ID line ends in a
#                      suffix bwa CANNOT handle (.1 .2 :a :b -- NOT /1,/2, which bwa
#                      auto-strips), else exit 1. Prints nothing; it is a predicate.
#                      Full scan, no early exit: absence of a bad suffix can only be
#                      proven by reading the whole file (a single file may mix
#                      conventions), and an early awk exit would SIGPIPE the upstream
#                      zcat and, under `set -o pipefail`, invert the gate.
#
set -euo pipefail

mode="${1:-strip}"

case "$mode" in
    strip)
        exec awk 'NR % 4 == 1 { sub(/[/.:][12ab]$/, "") } 1'
        ;;
    needs-normalizing)
        exec awk 'NR % 4 == 1 && /[.:][12ab]$/ { found = 1 } END { exit !found }'
        ;;
    *)
        echo "usage: $(basename "$0") [strip|needs-normalizing]" >&2
        exit 2
        ;;
esac
```

Note: `exec awk` replaces the shell, so awk's exit code becomes the script's exit
code directly (`set -e` no longer applies after `exec`). The `strip` awk program and
default mode are byte-for-byte the pre-existing behavior; only the dispatch and the
new mode are added.

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
docker run --rm --pull always -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest \
  bash testing/t/normalizeReadNames.t.sh
```
Expected: `PASS`

- [ ] **Step 5: Commit**

```bash
git add bin/normalizeReadNames.sh testing/t/normalizeReadNames.t.sh
git commit -m "Add needs-normalizing predicate mode to normalizeReadNames.sh

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: Gate normalization behind detection in `bwaMem`

**Files:**
- Modify: `modules/alignment.nf` (the `bwaMem` process `script:` block, currently lines ~40-62)

- [ ] **Step 1: Replace the `bwaMem` `script:` block**

Replace the entire `script:` triple-quoted block of the `bwaMem` process with the
following. (Leave the `container`, `input:`, `output:`, and `stub:` sections
untouched.) This preserves the existing line-continuation style: bwa argument lines
end with `\\` (literal backslash to the shell); the `samtools` chain lines end with
` && \` (Groovy line continuation), exactly as in the current file.

```groovy
  script:
      """
      set -euo pipefail

      if [ "$isPaired" = true ]; then
          if zcat -f $sample_1p | normalizeReadNames.sh needs-normalizing \\
             || zcat -f $sample_2p | normalizeReadNames.sh needs-normalizing; then
              # Read names carry a mate suffix bwa-mem2 cannot auto-strip
              # (.1/.2 from SRA, :a/:b from legacy Illumina; bwa only handles /1,/2).
              # Strip them on the fly so mem_sam_pe does not abort with
              # "paired reads have different names". Stripping forces bwa to pair
              # positionally, so first guard that the mate files are in lockstep;
              # fail loud rather than silently mispair.
              c1=\$(zcat -f $sample_1p | wc -l)
              c2=\$(zcat -f $sample_2p | wc -l)
              if [ "\$c1" -ne "\$c2" ]; then
                  echo "ERROR: mate files for ${sampleName} differ in length (\$c1 vs \$c2 lines); cannot positionally pair reads." >&2
                  exit 1
              fi
              bwa-mem2 mem \\
                  -t $params.bwaThreads \\
                  -R '@RG\\tID:${sampleName}\\tSM:${sampleName}\\tPL:ILLUMINA' \\
                  genomeIndex \\
                  <(zcat -f $sample_1p | normalizeReadNames.sh) \\
                  <(zcat -f $sample_2p | normalizeReadNames.sh) \\
                  | samtools collate -o output.bam - && \
                  samtools fixmate -m output.bam fix.bam && \
                  samtools sort -o sort.bam fix.bam && \
                  samtools markdup -r sort.bam result_sorted.bam
          else
              # Read names already conform (identical, or standard /1,/2 that bwa
              # auto-strips); hand the raw files straight to bwa untouched.
              bwa-mem2 mem \\
                  -t $params.bwaThreads \\
                  -R '@RG\\tID:${sampleName}\\tSM:${sampleName}\\tPL:ILLUMINA' \\
                  genomeIndex \\
                  $sample_1p \\
                  $sample_2p \\
                  | samtools collate -o output.bam - && \
                  samtools fixmate -m output.bam fix.bam && \
                  samtools sort -o sort.bam fix.bam && \
                  samtools markdup -r sort.bam result_sorted.bam
          fi
      else
          bwa-mem2 mem \\
              -t $params.bwaThreads \\
              -R '@RG\\tID:${sampleName}\\tSM:${sampleName}\\tPL:ILLUMINA' \\
              genomeIndex \\
              $sample_1p \\
              | samtools sort -o result_sorted.bam -
      fi
      """
```

- [ ] **Step 2: Verify Nextflow parses the module**

Run (stub mode compiles every process definition without needing real data/Docker):
```bash
nextflow run main.nf -entry processSingleExperiment -profile processSingleExperiment -stub-run
```
Expected: the DSL compiles and the run reaches process execution (stubs `touch`
outputs). A compile error in `bwaMem` would abort immediately with a Groovy parse
error naming `modules/alignment.nf`; that must NOT happen. (If the profile requires
inputs the stub run cannot satisfy and it fails *after* successful compilation, that
is acceptable for this step — the goal is to confirm the module parses. If Nextflow
is not installed locally, skip to Step 3 and rely on the reviewer's container run.)

- [ ] **Step 3: Verify the gate logic by extracting the paired branch**

Sanity-check the branch shape in isolation (independent of Nextflow variable
interpolation) using the real predicate script:
```bash
docker run --rm -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest bash -c '
  set -euo pipefail
  norm() { bin/normalizeReadNames.sh "$@"; }
  mk() { printf "%s\nACGT\n+\nFFFF\n" "$@"; }
  # conforming -> gate false -> raw path
  mk "@r/1" > c1.fq; mk "@r/2" > c2.fq
  if zcat -f c1.fq | norm needs-normalizing || zcat -f c2.fq | norm needs-normalizing;
     then echo "conforming: WRONG (took normalize path)"; exit 1; else echo "conforming: raw path OK"; fi
  # non-conforming -> gate true -> normalize path
  mk "@r.1" > n1.fq; mk "@r.2" > n2.fq
  if zcat -f n1.fq | norm needs-normalizing || zcat -f n2.fq | norm needs-normalizing;
     then echo "non-conforming: normalize path OK"; else echo "non-conforming: WRONG (took raw path)"; exit 1; fi
'
```
Expected:
```
conforming: raw path OK
non-conforming: normalize path OK
```

- [ ] **Step 4: Commit**

```bash
git add modules/alignment.nf
git commit -m "Gate read-name normalization on non-conformance in bwaMem

Detect mate suffixes bwa cannot handle before normalizing; pass
already-conforming paired samples to bwa untouched, and only assert
mate-length lockstep on the normalize path.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: Run the full test suite

**Files:** none (verification only)

- [ ] **Step 1: Run all repo tests in the container**

Run:
```bash
docker run --rm --pull always -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest bash -c '
  set -e
  for t in testing/t/*.jl; do julia "$t"; done
  python3 -m pytest testing/t/ -q
  bash testing/t/mergeBoth.t.sh
  bash testing/t/normalizeReadNames.t.sh
'
```
Expected: all Julia suites pass, pytest reports passed (e2e `test_mergeExperiments_e2e.py`
skips without `--run-dir`), `mergeBoth.t.sh` prints `PASS`, `normalizeReadNames.t.sh`
prints `PASS`.

- [ ] **Step 2: No commit** — verification only. If anything fails, return to the
  relevant task; do not paper over failures.

---

## Acceptance / follow-up (not automated here)

Matches PR #25's verification posture: the container smoke test on the exact failing
`IQ07` input remains the end-to-end proxy. A full cluster rerun is out of scope for
this plan. When updating the PR, note that conforming paired samples now bypass the
normalizer and the lockstep check entirely.
