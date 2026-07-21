# Conditional read-name normalization — design

**Date:** 2026-07-21
**Branch:** `fix-mate-readname-mismatch` (extends PR #25)
**Status:** approved, ready for implementation plan

## Problem

PR #25 fixes `bwa-mem2 mem` aborting paired alignment with
`[mem_sam_pe*] paired reads have different names` by feeding both mates through
`bin/normalizeReadNames.sh` (strips trailing `/1 /2`, `.1 .2`, `:a :b` from
read-ID lines) inside a process substitution, guarded by a mate-length
(lockstep) assertion.

As shipped, that machinery runs on **every** paired sample — including the
majority that already conform to the standard `/1`,`/2` suffix (or already have
identical mate names), which `bwa-mem2` handles natively with no help. The goal:
only pay for normalization on samples that actually need it, and leave
conforming samples' bytes untouched on their way to bwa.

## Core rule

`bwa-mem2 mem` pairs mate1/mate2 by **identical read name**, auto-stripping only
a trailing `/1` or `/2`. Therefore a paired sample is **already conforming** iff
no read-ID carries a suffix bwa cannot handle — i.e. none end in `.1 .2 :a :b`
(the `[.:][12ab]` set). A sample using only `/1`,`/2` (or no suffix) is
conforming and must be passed to bwa **raw**.

The detection set deliberately **excludes `/`**: `/1`,`/2` are bwa's job, so a
sample using the standard suffix is treated as conforming and skipped.

## Robustness constraint: mixed-convention files

PR #25's own motivating case (PlasmoDB `IQ07`) has a **single mate file that
mixes conventions** — a `090126:3:*` lane (identical-name, aligns fine)
interleaved with `90206:7/8:*` lanes (`:a`/`:b`, aborts bwa). Consequences:

1. Detection **cannot peek at only the first read** — a conforming first read
   can precede non-conforming reads later in the file. Detection must full-scan.
2. Once any non-conforming read is found, the **whole sample** is normalized
   (normalizing the already-conforming reads in the file is a harmless no-op for
   them — they have no `[/.:][12ab]$` suffix to strip).

## Design

### 1. `bin/normalizeReadNames.sh` — add a detection mode

One script, one source of truth for the suffix set, two modes:

- **`strip`** (default, unchanged behavior): reads FASTQ on stdin, strips
  `[/.:][12ab]$` from every 4th line, writes to stdout. Pure stream filter.
- **`needs-normalizing`**: reads FASTQ on stdin, exits **0** if any read-ID line
  (every 4th line) ends in `[.:][12ab]` (bwa-unhandled suffixes; excludes `/`),
  exits **1** otherwise. Prints nothing — it is a predicate.

**Full scan, no early `exit`.** Two reasons:
- Robustness — the absence of a bad suffix can only be proven by reading the
  whole file (mixed-convention files); there is no safe early-out on the
  conforming side.
- Correctness under `set -o pipefail` — an early `awk` `exit` inside
  `zcat -f … | normalizeReadNames.sh needs-normalizing` makes `zcat` die on
  SIGPIPE (exit 141), and `pipefail` would surface that as the pipeline status,
  **inverting the gate**. A full scan lets `zcat` finish cleanly (exit 0) so the
  `if` reads the predicate's real exit code.

Exit-code contract (`needs-normalizing`): `0` = needs normalizing, `1` =
conforming, `2` = usage error (unknown mode).

### 2. `modules/alignment.nf` — paired branch detects, then branches

```
if [ "$isPaired" = true ]; then
    if zcat -f $sample_1p | normalizeReadNames.sh needs-normalizing \
       || zcat -f $sample_2p | normalizeReadNames.sh needs-normalizing; then
        # non-conforming: positional pairing is now assumed, so guard lockstep
        c1=$(zcat -f $sample_1p | wc -l)
        c2=$(zcat -f $sample_2p | wc -l)
        if [ "$c1" -ne "$c2" ]; then
            echo "ERROR: mate files for <sample> differ in length ($c1 vs $c2 lines); cannot positionally pair reads." >&2
            exit 1
        fi
        bwa-mem2 mem … genomeIndex \
            <(zcat -f $sample_1p | normalizeReadNames.sh) \
            <(zcat -f $sample_2p | normalizeReadNames.sh) | samtools collate …
    else
        # conforming: hand raw files straight to bwa, untouched
        bwa-mem2 mem … genomeIndex $sample_1p $sample_2p | samtools collate …
    fi
else
    # single-end branch unchanged
fi
```

The two paired sub-branches share the same downstream
`samtools collate | fixmate | sort | markdup` pipeline; only the bwa input
differs (normalized process substitutions vs. raw file paths).

### 3. Lockstep assertion moves *inside* the normalize branch

The mate-length equality check exists **only because** stripping names to
identical forces bwa to pair positionally. A conforming sample pairs by real
read name and needs no length guard. Moving the assertion into the normalize
branch is strictly more correct than PR #25, which asserts lockstep even on
samples it never touches.

## Testing

Add `testing/t/normalizeReadNames.t.sh` (bash, following `mergeBoth.t.sh`
conventions: `mktemp -d`, `trap … EXIT`, `PASS`/`FAIL … exit 1`), run in the
`veupathdb/dnaseqanalysis:latest` container.

`strip` mode:
- strips `/1 /2`, `.1 .2`, `:a :b` from read-ID lines;
- leaves `keep/1extra` (suffix not at end) and `plain.name` untouched;
- only touches every 4th line (sequence/`+`/quality lines pass through).

`needs-normalizing` mode:
- conforming `/1`,`/2` file → exit **1**;
- identical-name file (no suffix) → exit **1**;
- `.1`/`.2` file → exit **0**;
- `:a`/`:b` file → exit **0**;
- mixed file (conforming reads + one `:a`/`:b` pair) → exit **0**;
- unknown mode → exit **2**.

## Cost trade-off

- **Conforming samples (common):** one full decompress-scan per mate for
  detection, then raw bytes to bwa. Replaces PR #25's two `wc -l` passes + inline
  awk — roughly a wash, with the payoff that conforming data is untouched and the
  spurious lockstep check is gone.
- **Non-conforming samples (rare):** detection scan + lockstep counts + the
  normalize pass — one extra pass vs. PR #25. Acceptable given rarity.

## Out of scope

- Tightening the strip regex itself (e.g. `:1`/`:2` false-positives on legacy
  tile coordinates) — pre-existing behavior of PR #25's normalizer, unchanged
  here; detection stays consistent with what `strip` acts on.
- End-to-end cluster rerun — same verification posture as PR #25 (container
  smoke test on the failing input as proxy).
