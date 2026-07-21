# Multiallelic locus reporting fixes — design

**Date:** 2026-07-21
**Branch:** `fix-multiallelic-locus-reporting`
**Component:** `bin/processSequenceVariations.jl`

## Motivation

Investigating `Pf3D7_01_v3:481838` (a haploid *P. falciparum* run,
`~/dnaseq_test/failed`) surfaced two independent defects at a multiallelic
locus. The input `merged.vcf.gz` there is a single multiallelic record
`G → T,C` with, among 216 samples: 65 `T`, 2 `C`, 59 covered reference
(synthesized from no-calls), 90 uncovered no-calls.

Two consumers of that locus disagree with the truth in two different ways:

1. **Genome browser (reads the split output VCF)** reports 61 and 124
   "homozygous for reference" in the `G→T` and `G→C` rows respectively — both
   inflated and mutually contradictory (true reference = 59).
2. **`variationFeature.dat`** reports `snp_major=T (0.5118)`,
   `snp_minor=C (0.0157)` — omitting the reference `G (0.4724)`, which is
   actually the *second-most-common* allele at the locus.

The two defects are in separate code paths and are fixed independently.

---

## Bug 1 — output VCF split fabricates reference for other-alt samples

### Root cause

`remap_sample_for_split` (the per-ALT biallelic splitter used only by
`write_vcf_entry`) remaps genotype allele indices as:

```julia
remap = idx -> idx == 0 ? 0 : (idx == target_alt_i ? 1 : 0)   # BUG: "other alt → 0"
```

When `G→T,C` is split into a `G→T` row and a `G→C` row, a sample genotyped as
`C` (index 2) in the `G→T` row is remapped `2 → 0` — i.e. **written as
homozygous reference** although it carries `C`, not the reference. Symmetrically
`T` samples are written as reference in the `G→C` row. This reproduces the
browser's 61 / 124 exactly (61 = 59 true-ref + 2 C-samples; 124 = 59 + 65
T-samples).

This is the same class of defect as bug (B) in the frequency-contract memory
(`bcftools norm --multi-overlaps 0` fabricating reference on split), which was
fixed at the *merge* step via `--multi-overlaps .`. This output writer
re-introduces it, hardcoded to `0`.

### Fix

A sample carrying a *different* ALT is **not** reference in this biallelic view;
it is **missing** (`.`), its real call living in the sibling split row (the
`--multi-overlaps .` convention). Change the remap so other alt indices → `.`:

```julia
remap = idx -> idx == 0 ? "0" : (idx == target_alt_i ? "1" : ".")
result[fi] = replace(gt, r"\d+" => m -> remap(parse(Int, m)))
```

`.`-slots and separators (`/`, `|`) are already preserved verbatim, so this is
ploidy-agnostic: a diploid compound het `1/2` splits to `1/.` and `./1` (no
fabricated reference), and a hom-alt `2/2` splits to `./.` in the `T` row.

### Expected result at 481838

| Row | variant | homozygous-ref | no-call |
|---|---|---|---|
| G→T | 65 | 59 | 92 (2 C + 90 uncovered) |
| G→C | 2  | 59 | 155 (65 T + 90 uncovered) |

### Explicitly out of scope

`AD`/`AO`/`RO`/`DP` in a split row still carry the full multiallelic vectors
(a `T` sample in the `G→C` row still shows its `T` reads in `AD`). Because that
sample is now `GT=.`, a per-record consumer that keys on `GT` ignores it.
Recomputing per-record `AD` (what full `bcftools norm` does) is a larger change
and is **not** part of this work. Confirmed the genome browser reads `GT`, so
the `GT` fix is sufficient for the reported symptom.

---

## Bug 2 — `variationFeature.dat` major/minor ignores the reference allele

### Root cause

`write_snp_feature` builds its major/minor ranking pools (`snp_keys`,
`indel_keys`) from **ALT alleles only**. The reference weight is siphoned into a
separate `ref_weight` accumulator (for `ref_allele_frequency`) and never
competes for a major/minor slot. Consequently `snp_major_allele` is the most
common *alt*, even when the reference is more common. In this run that misreports
**523,564 of 532,919 SNV rows (98.2%)**, where the reference is the true
majority allele but an alt is labelled `snp_major_allele`.

### Fix — per-class candidate pool includes the reference

For **each** variant class (SNP, indel), the ranking pool is that class's ALT
keys **plus** the reference key `(r, r)` whose reference string `r` matches the
class's alleles:

- SNP reference = the reference key whose ref-string equals the SNP alts'
  ref-string (e.g. `G`). Matching by **ref-string**, not length, so an MNP
  (equal-length multi-base substitution) keeps its multi-base reference in the
  SNP class.
- Indel reference = the reference key whose ref-string is the indel ref span
  (e.g. `ACA`), representing "no indel here".

Rank the pool by the existing key `(-weight, allele, ref)`. `major` = pool[1],
`minor` = pool[2]. Top-2 only: a third allele (e.g. `C` at 481838) drops off the
summary row but is still fully represented in `allele.dat`.

**The reference joins a class's pool only when that class has at least one ALT
allele.** A class with no alts stays empty (both major and minor blank) — a
locus with only the reference allele of a class produces no major/minor for that
class, exactly as today. This also keeps `variant_type` keyed on alt presence.

**Tie-break:** the existing sort key `(-weight, allele_string, ref_span)` is
applied to the combined pool. When the reference and an alt have equal weight,
the one whose allele string sorts first wins the higher slot (deterministic,
matches current behavior). Example: at a MIXED locus with SNP ref `A` (weight 1)
tied with SNP alt `G` (weight 1), `snp_major = A` (reference), `snp_minor = G`.

### Reference-wins-a-slot handling

- `snp_ref_allele` / `indel_ref_allele` columns are unchanged — always the
  genomic reference string.
- The slot (major or minor) holding the reference emits **empty**
  `genomic_hgvs` (no `g.481838G>G`). Alt slots keep normal `g.POS ref>alt`.
- The reference slot's `frequency` and `strain_count` come from the reference
  key's own `AlleleStat`. The reference key's `strains` set includes the
  synthetic reference strain (3D7), so at 481838 `snp_minor = G, freq 0.4724,
  count 60` — identical to `ref_allele_frequency` by construction. (This is +1
  versus a pure sample count; accepted for consistency with the existing
  `ref_allele_frequency` column and the frequency contract.)

### Unchanged

`ref_allele_frequency` (aggregate reference over the whole position),
`distinct_strain_count`, `called_strain_count`, `no_call_strain_count`,
`total_ploidy_count`, `het_strain_count`, and `variant_type` (still keyed on the
presence of SNP/indel **alt** alleles, computed before the reference is added to
the pools). Diploid ploidy-weighting is already correct: the reference key
accrues 2 weight units per diploid reference chromosome via
`chromosome_alleles`, so a diploid reference competes with its true 2N count.

### Expected result at 481838

`snp_ref_allele=G`, `snp_major=T (0.5118, 65)`, `snp_minor=G (0.4724, 60)`,
`snp_major_genomic_hgvs=Pf3D7_01_v3:g.481838G>T`, `snp_minor_genomic_hgvs=`
(blank). Across the run, the 98.2% ref-majority rows flip to `major=ref`,
`minor=the alt`, populating `minor` where it was previously empty.

---

## Testing

TDD, failing test first for each bug:

- `testing/t/handleVariantRecord.jl` — add cases: (a) split remap marks
  other-alt as `.` for haploid and diploid/compound-het; (b) `write_snp_feature`
  selects the reference as major/minor when it outranks alts, with blank HGVS on
  the reference slot; (c) multiallelic top-2 drops the third allele.
- `testing/t/test_mergeExperiments_e2e.py` — update existing assertions on
  major/minor columns and per-row genotype counts to the new semantics.

## Downstream / contract impact

Amend memory `project_variant_frequency_contract` and note that the web app must
now rank the reference into major/minor and expect a blank HGVS on a
reference-major/minor slot. Regenerate DAT expectations for the e2e suite.
