# Separate SNP and indel handling in mergeExperiments (via `-m both`)

**Date:** 2026-07-06
**Status:** Design — approved shape, pending spec review
**Branch target:** feature branch off `reduce-merge-outputs-snpeff-source`

## Problem

The website's variant record is locus-based (SNP + indel together), which is why
`mergeExperiments` merges per-experiment VCFs. But the merge step corrupts variant
classification:

`mergeVcfs` (modules/mergeExperiments.nf:44) runs `bcftools merge --merge all`,
which collapses SNP and indel alleles at the same start position into a **single
mixed multiallelic record** (e.g. `REF=A ALT=G,AT`). The annotator then has to
*guess* each row's class:

- `is_snp_record` (processSequenceVariations.jl:732) requires *all* alts to match
  REF length, so a mixed row is silently classified as "not a SNP record."
- `pick_snp_record` (line 2009) picks one record per locus and the output/CANN loop
  iterates only that record's alts, reading every strain's GT from it.

Consequences:

1. **Indels get no functional annotation.** They never reach `output.vcf`, so
   snpeff never sees them; they get no CANN entry in the mixed-row path.
2. **`downstream_of_frameshift` dead-ends for indels.** The flag is computed for
   indel variations (line 2210) but is never surfaced in CANN — an oversight.
3. **SNP/indel confusion** at mixed rows: GT for an indel variation is read from a
   record that may be the SNP record, and vice versa.

Empirically confirmed (bcftools 1.19):

```
-m all  (current):   chr1 100 A G,AT   1/1 2/2      # SNP + indel fused
-m both (proposed):  chr1 100 A G      1/1 ./.       # SNP record
                     chr1 100 A AT     ./. 1/1       # indel record, separate
```

## Key insight

The mixing is **not** inherent to merging atomized inputs — it is entirely an
artifact of the `-m all` flag. Per-experiment inputs already arrive atomized with
SNPs and indels on separate rows. Switching to `bcftools merge -m both` keeps
same-class alleles together and **never** merges across classes, producing
class-pure rows in a single merged VCF. No physical split into two files is needed.

## Invariant: at most 2 rows per locus

`-m both` packs *all* same-class alleles at a start position into one multiallelic
record each. Confirmed with 5 samples (2 different-length insertions + 1 deletion +
2 SNPs):

```
chr1 100 A    G,C            # ALL SNP alts, one record
chr1 100 ATG  ATTG,ATTTG,A   # ALL indels (padded to common REF), one record
```

So a locus yields **at most one multiallelic SNP record + one multiallelic indel
record** — never more, regardless of sample count or indel-length diversity. The
Julia loop's per-record iteration is therefore bounded at 2 records, each with N
alts. `record.alts` iteration already handles the N-alts dimension.

### REF padding caveat

When indels of different lengths co-occur, bcftools reconciles them to a **common
padded REF** (the longest span) and rewrites sibling ALTs accordingly
(`A>AT` becomes `ATG>ATTG` alongside a `ATG>A` deletion). Frame math is safe:
`len_diff` (alt_len − ref_len) is invariant under this padding, so
`build_cann_string`'s `len_diff % 3` classification still holds. But the indel
record the annotator sees has a multi-base REF with mixed insertion/deletion
siblings and `./.` for every strain not carrying a given allele. This is the shape
tests must exercise — not a clean single insertion. It also brushes against the
separately-deferred `matches_reference` overlapping-ref issue (out of scope here).

## Design

### Data flow

```
per-exp atomized VCFs ── bcftools merge -m both ── merged.vcf   (SNP & indel = separate rows)
                                                        │
              processSequenceVariations.jl (rows are class-pure; class is explicit)
                    ├─ output.vcf  (SNP rows + indel rows, separate)
                    │       └─ snpEff (one run) ─ parseSnpEffAnnotations ─ snpeff.dat
                    └─ allele / variation / product .dat   (locus-joined in the single pass, as today)
```

### Changes

| # | Change | Location | Size |
|---|--------|----------|------|
| 1 | `--merge all` → `--merge both` | modules/mergeExperiments.nf:44 | trivial |
| 2 | Emit indel records in output/CANN path; iterate **both** class records per locus, reading each variation's GT from its own class-record | `handle_variant_record!` output loop + `collect_cann_entries_for_annotation` | **bulk** |
| 3 | Surface `downstream_of_frameshift` for indels in CANN as a compound effect | `build_cann_string` (~line 1943) | small |
| 4 | Ensure effect terms map (incl. `downstream_frameshift`) | parseSnpEffAnnotations.py `CANN_EFFECT_MAP` | small |
| 5 | Spike + tests for `./.` attribution and padded multiallelic indel record | testing/t | moderate |
| 6 | Update golden `.dat` / snpeff expectations | testing/t | moderate |

**Not touched:** the `snpEff` process runs once and is unchanged (snpeff handles
indels natively). `parseSnpEffAnnotations.py` already carries the indel effect
vocabulary (`frameshift`, `inframe_insertion_unnormalized`,
`inframe_deletion_unnormalized`) and already emits multiple rows per locus. The
`.dat` locus-grouping logic, the frameshift SQLite DB, and per-experiment
`filterAndSplitVcf` are unchanged. No new Nextflow processes; no file split; no
`.dat` re-merge step.

### Output-loop change (the bulk — #2)

Today: `record = pick_snp_record(records)`; the output loop iterates that single
record's alts and reads all GTs from it. With `-m both` a locus has up to two
records (SNP + indel), each with its own sample columns. The loop must iterate
**both** records and attribute each variation's GT from its own class-record. This
is the same code region as the original bug, which is why it is de-risked with a
spike first.

`build_cann_string` already computes indel frame effects
(frameshift / inframe_insertion / inframe_deletion, lines 1943–1952) — currently
dead code for the mixed-row path. This design finally exercises it.

### Indel + downstream-of-frameshift representation (#3)

Emit both terms as a compound effect, e.g. `frameshift&downstream_frameshift` or
`inframe_deletion&downstream_frameshift`. The structural term states what the indel
does to length; `downstream_frameshift` carries the context that the reading frame
is already broken upstream. Compound effects already use `&` and
`parseSnpEffAnnotations.py` already splits on it into one row per effect, so no
parser structural change is needed beyond mapping the term.

## `./.` GT semantics (primary risk)

`-m both` changes GT patterns: a strain carrying the SNP shows `./.` on the indel
record and vice versa. The allele-attribution code (`build_variations_from_record`,
`fill_missing_coverage_gt`) must treat `./.` on a class-record as **"this strain has
no variant of this class here"** — it must NOT fabricate a reference call nor
trigger coverage-fill that invents data.

## Testing

**Spike first (before broad implementation):** a focused Julia test that builds a
single-locus fixture with multiple strains, ≥2 different-length insertions **and** a
deletion (so the indel record is multiallelic with a padded REF), plus at least one
SNP. Assert:

- allele attribution is correct with `./.` rows (no fabricated ref calls / bogus
  coverage fill);
- indel CANN frame effect is correct on the padded multiallelic record;
- `downstream_of_frameshift` surfaces on the indel CANN when applicable.

**Then:** unit tests for indel CANN incl. the compound `downstream_frameshift`
term (`testing/t/handleVariantRecord.jl`); python parser tests for the new/mapped
terms (`testing/t/test_parseSnpEffAnnotations.py`); refresh golden `.dat` and
snpeff expectations; run the Nextflow test profile.

## Out of scope

- The deferred `matches_reference` overlapping-ref bug.
- Any change to per-experiment `processSingleExperiment` VCF processing.
- SnpEff's own indel normalization vs. Julia's (already handled via the
  `_unnormalized` effect terms).
