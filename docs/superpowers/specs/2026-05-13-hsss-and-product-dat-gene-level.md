# HSSS Binary Strain Files + Output File Restructure

**Date:** 2026-05-13  
**Status:** Approved

## Overview

Three related changes to `processSequenceVariations.jl` and `modules/mergeExperiments.nf`:

1. **variationFeature.dat** — strip to purely genomic columns; one row per position; drop `transcript_id` and all CDS-specific columns
2. **transcript_product.dat** — new file replacing both `product.dat` and `cache.tsv`; per (position, transcript, codon); also serves as cache input for subsequent runs
3. **HSSS binary strain files** — new output for HighSpeedSnpSearch; 4 read-frequency cutoffs

---

## Part 1: variationFeature.dat — Purely Genomic

### Problem

`write_snp_feature` is called once per transcript inside the annotation loop, producing multiple rows per position and mixing genomic-level and transcript-level data in the same file.

### Solution

Move `write_snp_feature` outside the annotation loop — called once per position. Drop all CDS/transcript-specific columns; those move to `transcript_product.dat`.

### Dropped Columns

| Column | Reason |
|---|---|
| `transcript_id` | Per-transcript; caused multiple rows |
| `has_nonsynonymous_allele` | Derivable from transcript_product.dat |
| `major_product`, `minor_product` | Derivable from transcript_product.dat |
| `has_stop_codon` | Derivable from transcript_product.dat |
| `ref_codon` | Per-transcript; moves to transcript_product.dat |
| `pos_in_cds` | Per-transcript; moves to transcript_product.dat |
| `downstream_of_frameshift_strain_ids` | Per-transcript; moves to transcript_product.dat |

### Resulting Columns (14)

```
location  seq_id  reference_strain  ref_allele
major_allele  minor_allele  major_allele_strain_count  minor_allele_strain_count
major_allele_frequency  minor_allele_frequency
distinct_strain_count  distinct_allele_count  total_ploidy_count  is_coding
```

`is_coding` is position-level: `1` if any overlapping transcript is coding, else `0`.
Allele counts/frequencies are computed from variations, not from annotations — no transcript dependency.

---

## Part 2: transcript_product.dat

### Purpose

Replaces both `product.dat` (per-transcript codon/AA rows) and `cache.tsv` (per-transcript CDS position cache). The first four columns preserve the existing cache format so `parse_cache_tsv_record` requires minimal changes.

The file header starts with `#` so the cache reader skips it correctly.

### When Written

Per (position, transcript) — inside the annotation loop, only for coding positions with variation (unchanged from current `write_product_file` call site).

### Columns (12)

```
#location  seq_id  transcript_id  pos_in_cds  pos_in_protein
codon  pos_in_codon  count  product  matches_ref_codon  matches_ref_product
downstream_of_frameshift_strain_ids
```

- `pos_in_protein` = `div(pos_in_cds - 1, 3) + 1` (1-based AA position; computed at write time)
- `downstream_of_frameshift_strain_ids` — same format as current variationFeature.dat: `{id1,id2,...}` or `""`; same value repeated for all codon rows of a given (position, transcript)
- One row per unique expanded codon per (position, transcript)

### As Cache Input

`parse_cache_tsv_record` reads columns 1–4 (`location, seq_id, transcript_id, pos_in_cds`) and ignores the rest — backwards compatible. The `--cache_file` CLI argument now accepts a previous run's `transcript_product.dat`.

The separate `--output_cache` CLI argument and `output_cache.tsv` output are eliminated.

### Nextflow Changes

| Before | After |
|---|---|
| `path 'product.dat', emit: productFile` | `path 'transcript_product.dat', emit: transcriptProductFile` |
| `path 'output_cache.tsv', emit: outputCache` | *(removed)* |
| `publishDir ... 'product.dat'` | `publishDir ... 'transcript_product.dat'` |
| `publishDir ... 'output_cache.tsv'` | *(removed)* |
| `--output_cache output_cache.tsv` in script | *(removed from CLI call)* |

---

## Part 3: HSSS Binary Strain Files

### Background

Replaces the Perl script `hsssCreateStrainFiles`. All logic is added directly to `processSequenceVariations.jl` — no intermediate TSV needed since all variation data is already in memory.

### Output Structure

Four directories: `hsss_readFreq20/`, `hsss_readFreq40/`, `hsss_readFreq60/`, `hsss_readFreq80/`

Each contains:

| File | Format | Description |
|---|---|---|
| `strainIdToName.dat` | TSV `index\tname` | Strain integer index → name |
| `contigIdToSourceId.dat` | TSV `seq_index\tsource_id` | Sequence integer index → source_id |
| `referenceGenome.dat` | Binary | Reference strain records |
| `{index}` (one per non-ref strain) | Binary | Per-strain records |

Reference strain = index `1`. Non-reference strains = indices `2…N` in `ctx.all_strains` order.

### Binary Record Format (8 bytes, matches Perl `pack("slcc", ...)`)

| Field | Julia type | Bytes | Description |
|---|---|---|---|
| `seq_index` | `Int16` | 2 | Sequence integer index |
| `location` | `Int32` | 4 | 1-based genomic position |
| `allele_code` | `Int8` | 1 | 0=unknown, 1=A, 2=C, 3=G, 4=T |
| `product_code` | `Int8` | 1 | ASCII of AA char, or 0 if unknown |

### Product Encoding

Across **all** annotations at a position (all transcripts, all genes):
- `unique(all_products)` has exactly one element → `Int8(codepoint(only(unique_products)))`
- Any disagreement (including overlapping opposite-strand genes) → `Int8(0)`

HSSS has a single product slot per record and cannot represent the overlapping-gene case.

### Allele Encoding

```julia
const HSSS_ALLELE_CODE = Dict('A'=>1,'a'=>1,'C'=>2,'c'=>2,'G'=>3,'g'=>3,'T'=>4,'t'=>4)
```

Any other base (IUPAC ambiguity, deletion, etc.) → `0`.

### Single-Pass Strategy

Open all 4 × (N_strains + 2) file handles before the main position loop.

For each position, for each of the 4 cutoffs:
- Strain "passes" if `parse(Float64, v.percent) >= cutoff`
- If no non-reference strain passes and differs from reference → skip this cutoff/position
- Reference strain → write to `referenceGenome.dat`
- Non-reference strain that passes + has variation → write record (all alleles if het)
- Non-reference strain that passes + matches reference with same product → skip
- Non-reference strain absent from this position → write unknown `(seq_idx, location, 0, 0)`

`hsss_seq_index` increments when `seq_id` changes; appended to `contigIdToSourceId.dat`.

### Nextflow Additions to `processSeqVars`

```nextflow
publishDir "$params.outputDir", mode: "copy", pattern: 'hsss_readFreq20'
publishDir "$params.outputDir", mode: "copy", pattern: 'hsss_readFreq40'
publishDir "$params.outputDir", mode: "copy", pattern: 'hsss_readFreq60'
publishDir "$params.outputDir", mode: "copy", pattern: 'hsss_readFreq80'

output:
  path 'hsss_readFreq20', emit: hsssReadFreq20
  path 'hsss_readFreq40', emit: hsssReadFreq40
  path 'hsss_readFreq60', emit: hsssReadFreq60
  path 'hsss_readFreq80', emit: hsssReadFreq80

stub:
  mkdir hsss_readFreq20 hsss_readFreq40 hsss_readFreq60 hsss_readFreq80
```

---

## Implementation Order

1. Restructure `write_snp_feature` — move call outside annotation loop; drop CDS columns
2. Update variationFeature.dat header (14 columns)
3. Rename `write_product_file` → `write_transcript_product` — add `pos_in_cds`, `pos_in_protein`, `downstream_of_frameshift_strain_ids` columns
4. Merge `cache_fh` into `transcript_product_fh`; remove `open`/`close` of separate cache file
5. Update `parse_cache_tsv_record` — handle ≥4 columns (already does); update header comment
6. Remove `--output_cache` CLI argument; update `open_output_writers`
7. Rename output file references: `product.dat` → `transcript_product.dat`
8. Add `HSSS_ALLELE_CODE` constant and product/allele encoding helpers
9. Add HSSS file handle initialization (open files, write index maps)
10. Add HSSS write calls inside main position loop (all 4 cutoffs, single pass)
11. Add HSSS teardown
12. Update `modules/mergeExperiments.nf` — remove cache output, rename product output, add HSSS outputs
13. Update Julia unit tests — `write_snp_feature` (column count), `write_transcript_product` (new columns), HSSS binary
14. Update e2e tests — variationFeature.dat column count, transcript_product.dat columns, HSSS directories

---

## Tests

- **Julia unit:** `write_snp_feature` emits exactly 14 tab-separated columns, one row per position
- **Julia unit:** `write_transcript_product` emits correct `pos_in_protein`, `downstream_of_frameshift_strain_ids`
- **Julia unit:** cache reader correctly reconstructs annotations from transcript_product.dat (reads cols 1–4)
- **Julia unit:** HSSS binary — 2-strain, 3-position fixture; verify byte layout for each cutoff
- **e2e:** variationFeature.dat has 14 columns, one row per variant position
- **e2e:** transcript_product.dat has 12 columns; header starts with `#`
- **e2e:** All 4 `hsss_readFreqXX` directories present with expected files
