# Complex/Multiallelic Variant Frequency Over-Counting — Design

**Date:** 2026-07-07
**Status:** Design approved; ready for implementation plan
**Affects:** `modules/mergeExperiments.nf` (`mergeVcfs`), `bin/processSequenceVariations.jl`
**Outputs corrected:** `variationFeature.dat`, `allele.dat`

## Problem

At loci where a single strain's diploid genotype is spread across **multiple VCF
records**, the merge pipeline counts allele frequencies per-record instead of
per-chromosome, producing wrong `variationFeature.dat` / `allele.dat` numbers.
Two distinct root causes were found during QA of the `merge` test dataset
(`/home/jbrestel/dnaseq_test/merge/output`).

### Bug A — complex-variant decomposition inflates the denominator (locus LmjF.01:85879)

FreeBayes calls one complex allele `TGT>C` (`1/1`, `AO=6`, `DP=6`) — a SNP
`T>C` at 85879 **plus** a 2 bp deletion of `GT` at 85880–81 on one haplotype.
`bcftools norm -a` (in `filterAndSplitVcf`, `modules/snp.nf:89`) faithfully
decomposes it into two records:

```
85879  T   > C    1/1   (SNP primitive)
85879  TGT > T    1/1   (deletion primitive)
```

Both records are `1/1`, so the same diploid strain contributes its full ploidy
(2) to the SNP class **and** again to the indel class. The reference strains are
counted once. Result: the shared denominator inflates.

| field | before | correct |
|---|---|---|
| `total_ploidy_count` | 13 | 9 |
| `ref_allele_frequency` | 0.3846 | 0.5556 |
| snp `C` freq | 0.3077 | 0.4444 |
| indel `del` freq | 0.3077 | 0.4444 |

The per-allele numerators (`C_weight=4`, `del_weight=4`, `ref_weight=5`) are
**correct** — only the denominator is wrong.

### Bug B — `1/2` multiallelic split fabricates reference alleles (locus LmjF.01:8962)

FreeBayes calls Seidman751 as a genuine compound heterozygote:
`T>TA,TAA  GT=1/2  AD=0,7,5  RO=0` — one chromosome carries `+A`, the other
`+AA`, **zero reference chromosomes**. `mergeVcfs` runs `bcftools norm -m -any`,
which splits the `1/2` into two biallelic records and — because its default is
`--multi-overlaps 0` — rewrites the complementary allele as reference:

```
T>TA    GT=1/0    (fabricated reference on the second chromosome)
T>TAA   GT=0/1    (fabricated reference on the first chromosome)
```

Downstream then counts 2 fabricated reference copies for a strain that has none.

| field | before | correct | why wrong |
|---|---|---|---|
| `total_ploidy_count` | 11 | 9 | Seidman counted as ploidy 4 |
| `ref_allele_frequency` | 0.8182 | 0.7778 | +2 fabricated ref copies |
| `het_strain_count` | 2 | 1 | counts records, not strains |
| TA / TAA freq | 0.0909 | 0.1111 | inflated denominator |
| allele.dat ref strain set | `{1,2,3,4,5}` | `{1,2,4,5}` | **Seidman falsely listed as carrying reference** |

`allele.dat` literally asserts Seidman carries a reference allele at a locus
where it does not.

### The VCF is correct — no VCF change needed

Both bugs are **downstream accounting artifacts**, not VCF errors:

- The raw FreeBayes complex call `TGT>C` is one `1/1` haplotype (6 reads); its
  atomization into SNP + deletion is a faithful primitive decomposition.
- The raw `1/2` compound het is one call (7 + 5 reads, `RO=0`); the fabricated
  reference is introduced by `bcftools norm -m -any`'s default fill, not by
  FreeBayes.

The `1M2D10M1X7M` CIGAR in the atomized records is a stale FreeBayes INFO field
not recomputed by `bcftools norm -f`; it is unreliable and must be ignored — the
`REF/ALT` fields are the source of truth.

### Common thread

One strain's diploid genotype is spread across multiple records, and per-record
accounting (a) double-counts ploidy (Bug A) and (b) fabricates reference alleles
(Bug B). A single split record `0/1` is indistinguishable from a genuine ref/alt
het, so the merge representation must preserve the distinction; the counter must
then honor it.

## Design

Two coupled changes that **must land together** — each alone breaks the other
case (see Coupling below).

### Change 1 — merge representation (`modules/mergeExperiments.nf`, `mergeVcfs`)

Add `--multi-overlaps .` to the `bcftools norm -m -any` invocation
(`modules/mergeExperiments.nf:41`):

```
bcftools norm -m -any --multi-overlaps . -Oz -o "${vcf%.vcf.gz}.norm.vcf.gz"
```

A split `1/2` then becomes `1/.` + `./1` — each record speaks for exactly one
chromosome, with the other chromosome's allele marked missing (`.`) rather than
fabricated reference (`0`). Verified in the `veupathdb/dnaseqanalysis:1.1.0`
container (bcftools 1.19). No effect on Bug A (two separate records, never a
multiallelic split).

### Change 2 — Julia counter (`bin/processSequenceVariations.jl`)

**2a. GT parsing must honor missing slots.** Today `gt_to_base` (line 1014) and
`nonref_alt_alleles` (line 1049) return `""` whenever *either* GT slot is `.`,
and `build_variations_from_record` (line 1150) then skips the record entirely.
Under Change 1 that would silently drop the TA/TAA insertions (undercount). Fix:
a half-missing GT (`1/.`) counts its present allele(s) as chromosome copies and
treats each `.` slot as **missing** — contributing nothing to any allele and not
dropping the record.

**2b. `aggregate_locus_alleles` — count per non-missing chromosome slot; set
denominator to distinct-strain ploidy.**

- **Numerator:** for each variation, add one unit of weight per **non-missing**
  GT slot to `(reference_span, allele)`; a `0` slot adds to the reference key, an
  ALT slot adds to that ALT key, a `.` slot adds nothing. This removes the het
  branch's unconditional `add!(ref, ref, 1, …)` that fabricates reference for
  missing slots.
- **Denominator (`total`):** `Σ ploidy` over **distinct strains** (each strain's
  `gt_ploidy` counted once), replacing the current `total += v.ploidy`
  accumulation per record.

The single `total` returned by `aggregate_locus_alleles` feeds every dependent
value — `total_ploidy_count` (line 1567), `ref_allele_frequency` (1535), snp/indel
major+minor frequencies (1541/1546), and `allele.dat` per-allele frequency
(1598) — so both output files are corrected consistently.

**2c. `het_strain_count`** must count distinct strains with a het call, not
variation records (today: `count(v -> !isempty(v.alt_allele), variations)`), so a
split compound het counts as one strain.

### Coupling

| | Change 1 only | Change 2 only | Both |
|---|---|---|---|
| Bug A (85879) | unchanged (still 13) | fixed (9) | fixed |
| Bug B (8962) | records become `1/.`, current code **drops them** → insertions vanish | ref still fabricated → ref freq → 1.0 (worse) | fixed |

## Net effect

| locus | field | before | after |
|---|---|---|---|
| 85879 | total_ploidy / ref / C / del | 13 / .3846 / .3077 / .3077 | 9 / .5556 / .4444 / .4444 |
| 8962 | total_ploidy / ref / TA / TAA | 11 / .8182 / .0909 / .0909 | 9 / .7778 / .1111 / .1111 |
| 8962 | Seidman in ref strain set | yes (wrong) | no |
| 8962 | het_strain_count | 2 | 1 |

Simple single-record loci (one record per strain) are unaffected: each strain is
counted once either way.

## Frequency contract (for web-app reuse on sample subsets)

This is the canonical definition the web app must replicate when computing
frequencies over any subset of samples:

1. **Chromosome-copy model.** Each strain contributes `ploidy` chromosome copies
   at a covered locus (diploid = 2; the synthetic reference strain = 1). A
   strain is counted **once**, even when its genotype is spread across multiple
   co-located records (complex decomposition or split multiallelic).
2. **Denominator** `total_ploidy_count` = `Σ ploidy` over the distinct strains
   present (called or covered-reference) in the subset, including the synthetic
   reference strain.
3. **Numerator** for an allele = number of chromosome copies carrying it =
   count of non-missing GT slots equal to that allele across the strain's
   records. A `.` (missing) slot contributes to no allele.
4. **Reference frequency** = reference copies / denominator. A strain with a
   genuine compound het (`1/2`) contributes **0** reference copies.
5. **Per-class frequencies may sum to more than 1.0** at complex/MIXED loci,
   because one chromosome can be simultaneously a SNP and an indel. This is
   correct, not a normalization error. `reference + Σ SNP alts + Σ indel alts`
   equals 1.0 only at loci with no co-located complex variant.

## Blast radius

- Only multi-record-per-strain loci change: complex decompositions (Bug A) and
  split multiallelic hets (Bug B). Quantify exact counts across the test set
  during implementation (`merge/output` had 69 multi-record positions of 1504).
- `snpEff` and consensus/`CANN` steps consume the merged VCF; verify they
  tolerate `1/.` / `./1` genotypes (spot-check during e2e verification).

## Testing

- Julia unit tests (`testing/t/handleVariantRecord.jl`): add cases for
  (a) complex decomposition (same strain `1/1` in a SNP and an indel record),
  (b) compound het `1/2` split into `1/.` + `./1`,
  (c) half-missing GT parsing in `gt_to_base` / `nonref_alt_alleles`.
- Relax any existing assertion that `ref + snp + indel == 1.0` to hold only for
  non-complex loci.
- Update fixture expected values for `total_ploidy_count` / frequencies at
  affected MIXED loci.
- e2e: re-run `mergeExperiments` on the test dataset; assert the 85879 and 8962
  rows match the "after" tables above.

## Out of scope

- No change to FreeBayes calling, `filterAndSplitVcf` atomization, or the
  published per-sample VCF.
- No change to the per-class `variationFeature.dat` column schema
  (see `2026-07-06-variationfeature-per-class-schema-design.md`).
