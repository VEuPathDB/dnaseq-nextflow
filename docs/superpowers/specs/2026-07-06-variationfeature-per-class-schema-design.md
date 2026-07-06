# variationFeature.dat per-class (SNP/indel) schema + allele identity fix

**Date:** 2026-07-06
**Status:** Implemented and merged to `reduce-merge-outputs-snpeff-source`. Unit tests green (`handleVariantRecord.jl`); e2e green (`test_mergeExperiments_e2e.py` — 96/96 against a fresh `mergeExperiments` run 2026-07-06). Both the deferred `matches_reference` overlapping-ref bug and Finding 2 (SNP+deletion allele collapse) are fixed and verified on real data at `LmjF.01:13850` (deletion `A` now a distinct `matches_reference=0` row with correct `delCA` HGVS, no longer collapsed into reference).
**Branch target:** feature branch off `reduce-merge-outputs-snpeff-source`

## Context

After `bcftools merge -m both` (branch `split-snp-indel-merge-both`, merged), a locus
can carry a SNP record and an indel record as separate VCF rows. But the `.dat`
writers still summarize a locus with a single shared allele model:
`write_snp_feature` / `write_allele_file` (`bin/processSequenceVariations.jl`) take
`ref_allele = variations[1].reference` (one arbitrary record's ref) and key alleles
by **bare base string**. At a mixed locus this conflates two biological meanings of
the same string — e.g. deletion product `A` (ref `ACA`) collides with SNP reference
`A` — producing the deferred `matches_reference` bug (overlapping multi-length records at
one POS) and Finding 2 (deletion allele collapses into the reference allele).

`variationFeature.dat` is the webapp's per-locus table: **one row = one locus = one
record page**, and its columns are what a user sees in a tabular list of loci.
`allele.dat` is QA-only. These files map directly to a (not-yet-created) postgres
schema, so we are free to choose columns.

**Goal:** represent SNP and indel information for a locus without conflation, fixing
both deferred bugs at their shared root, while keeping one row per locus.

## Decisions (from brainstorming)

1. **Wide flat `variationFeature`** — one row per locus: locus-invariant columns +
   `snp_*` family + `indel_*` family. Single-class loci leave the other family empty.
2. **Locus-wide shared frequency denominator** — one ploidy-weighted allele space per
   locus over all called strains; `reference` + all SNP alts + all indel alts sum to
   1.0. Per-class frequencies are read from this shared space so they are directly
   comparable on a record page.
3. **Allele identity = `(ref_span, allele)` tuple**, not the bare string — the actual
   bug fix.
4. **Reference counted once**, at locus level (`ref_allele_frequency`). Per-class
   `major`/`minor` are chosen among that class's **alt** alleles only.
5. **Major + minor per class.** Full per-allele detail remains in `allele.dat`.
6. **Derive flags, don't store** — no `major_differs_from_reference` / `is_singleton`
   columns; downstream derives them (`snp_major_allele != snp_ref_allele`;
   `strain_count == 1`).

## variationFeature.dat schema (31 columns)

**Locus-invariant (always populated):**
```
location, seq_id, reference_strain, is_coding, variant_type,
distinct_strain_count, called_strain_count, no_call_strain_count, call_rate,
total_ploidy_count, ref_allele_frequency, het_strain_count
```
- `variant_type ∈ {SNV, INDEL, MIXED}` (unchanged classification).
- `ref_allele_frequency` = ploidy-weighted fraction of the locus-wide space that is
  the reference genome. Well-defined even though the ref *string* differs per class.
- `het_strain_count` = number of strains with a het call at the locus
  (`count(v -> !isempty(v.alt_allele), variations)`), unchanged from today.

**SNP family — `snp_*` (empty when the locus has no SNP):**
```
snp_ref_allele, snp_major_allele, snp_major_allele_frequency, snp_major_allele_strain_count,
snp_minor_allele, snp_minor_allele_frequency, snp_minor_allele_strain_count,
snp_major_genomic_hgvs, snp_minor_genomic_hgvs
```

**Indel family — `indel_*` (empty when the locus has no indel):**
```
indel_ref_allele, indel_major_allele, indel_major_allele_frequency, indel_major_allele_strain_count,
indel_minor_allele, indel_minor_allele_frequency, indel_minor_allele_strain_count,
indel_major_genomic_hgvs, indel_minor_genomic_hgvs, indel_frame_effect
```
- `indel_frame_effect ∈ {frameshift, inframe_insertion, inframe_deletion}` for the
  major indel allele (length delta vs `indel_ref_allele`).

**Removed** (the collapsing columns): single `ref_allele`, `major_allele`,
`minor_allele`, `major_allele_*`, `minor_allele_*`, `major_differs_from_reference`,
`is_singleton`. Their meaning moves into the per-class families or is derived.
`het_strain_count` is retained as a locus-invariant column.

## Computation

1. **Build the locus-wide weighted allele map** from all variations (both records'
   variations, combined as today), ploidy-weighted (reuse `compute_allele_weight_map`
   semantics), but **keyed by `(reference, base)` tuples** rather than bare `base`.
   Denominator `total_ploidy_count` includes the synthetic reference strain, as today.
2. **Bucket each key:**
   - reference bucket: keys where `base == reference` → aggregated into
     `ref_allele_frequency`.
   - SNP alts: `base != reference && length(base) == length(reference)`.
   - indel alts: `length(base) != length(reference)`.
   (`v.alt_allele` het handling: a het `1/2` strain splits ploidy weight across its two
   alt keys via `compute_allele_weight_map`, landing in whichever class each alt belongs
   to — no special case.)
3. **Per class**, rank that class's alt keys by weight → `major` (top), `minor`
   (second). `*_strain_count` = distinct strains carrying that key. `*_frequency` =
   `weight / total_ploidy_count`. `snp_ref_allele` / `indel_ref_allele` = the
   `reference` field of that class's keys (one value per class under `-m both`, since a
   locus has ≤1 record per class).
4. **genomic HGVS** per allele via `allele_genomic_hgvs` using the allele's own ref
   span (from `build_allele_refs`, already keyed by allele — extend to key by
   `(ref,base)` consistently). `indel_frame_effect` from `length(major_indel) -
   length(indel_ref_allele)` mod 3.

## allele.dat (QA) — same fix

Re-key `allele_entries` / `allele_counts` by `(reference, base)` so the QA file stops
collapsing. `matches_reference` becomes per-row correct: `1` iff that row's
`base == reference`. Columns otherwise unchanged (still per (locus, allele) with its
own ref span + genomic_hgvs). This closes the deferred `matches_reference` bug.

## Components / writer refactor

`bin/processSequenceVariations.jl`:
- New helper: build the locus-wide `(ref,base)`-keyed weighted map + strain counts
  once, returned in a small struct, so `write_snp_feature` and `write_allele_file`
  share one correct aggregation (removes today's duplicated, divergent counting).
- `write_snp_feature` restructured: aggregate → partition into ref/snp/indel →
  emit the 30 columns. Header updated.
- `write_allele_file`: re-key by `(ref,base)`; `matches_reference` per row.
- `build_allele_refs` / `allele_genomic_hgvs`: key by `(ref,base)` consistently.

Keep each unit small and independently testable (the aggregation helper is pure and
unit-testable without writers).

## Testing

- **Unit** (`testing/t/handleVariantRecord.jl`, capturing writer `IOBuffer`s via
  `make_intergenic_ctx`): at a mixed SNP+deletion locus (the `LmjF.01:13850` shape:
  SNP `A>G`, deletion `ACA>A`, plus reference strains), assert:
  - `variationFeature`: `variant_type==MIXED`; `snp_*` and `indel_*` families both
    populated with correct ref spans (`A` vs `ACA`); the deletion is **not** under the
    SNP family and **not** mislabeled reference; per-class + `ref_allele_frequency` sum
    to 1.0; `indel_frame_effect` correct.
  - `allele.dat`: deletion `A` and reference are distinct rows; `matches_reference`
    correct on each (deletion row `0`, reference row `1`).
  - SNP-only and indel-only loci leave the opposite family empty.
- **e2e** (`testing/t/test_mergeExperiments_e2e.py --run-dir <run>`): update the
  `variationFeature`/`allele` column-count + column-name assertions to the new schema;
  add invariants (per-class frequency ranges, families empty iff class absent). Re-run
  against a real run before merge (needs Docker + a completed `mergeExperiments`).

## Out of scope

- Downstream GUS table + webapp query schema (built later from these files).
- `transcript_product.dat`, `snpeff.dat`, HSSS files (HSSS is SNP-only and already
  correct with separated rows).
- Any change to `-m both` / annotator record handling (already merged).
