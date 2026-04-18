# Het Indel Product Calling — Design Spec

**Date:** 2026-04-17  
**Status:** Approved

## Problem

Heterozygous indels cannot be represented in a single-sequence consensus FASTA because the two alleles differ in length. The current pipeline masks het indel positions with `N`, the same character used for zero-coverage positions. Homozygous indels are emitted as the called ALT allele, causing the consensus FASTA to diverge in length from the reference. This means:

1. Downstream CDS extraction must track coordinate shifts introduced by applied hom indels.
2. Het indel positions yield `N` → unknown amino acids, so no products are reported.
3. No products are reported for het indels, even though het SNPs report both products via IUPAC expansion.

## Goal

Report both REF-allele and ALT-allele products for heterozygous indels, consistent with how het SNPs are handled today. Simplify the coordinate system by keeping the consensus FASTA permanently reference-length.

## Core Principle

The consensus FASTA is a SNP + coverage artifact only. It never encodes indels.

| Position type | FASTA | DB |
|---|---|---|
| Low / no coverage | `N` | — |
| SNP | IUPAC code | — |
| Hom insertion | ref bases | yes |
| Hom deletion | ref bases | yes |
| Het insertion | ref bases | yes |
| Het deletion | ref bases | yes |

All indel information — zygosity, REF allele, ALT allele — lives exclusively in the `genomic_indels` table. `makeCodingData.jl` applies indels at product-calling time from the DB. The FASTA coordinate system is always reference space; no shifting required.

## Design

### Signal layer — `makeConsensusFastaFromGvcf.py`

Two changes to `build_consensus`:

1. **Hom indel branch** (currently lines 99–102): replace `segments.append(alleles[0])` with `segments.append(ref_seq[pos:pos + len(v.REF)])`. The REF span is preserved; the FASTA length never diverges from the reference.

2. **Het indel branch** (currently lines 104–107): replace `segments.append('N' * len(v.REF))` with `segments.append(ref_seq[pos:pos + len(v.REF)])`. `N` is reserved for low/no coverage only.

### Data layer — `genomic_indels` table

**`findValues.pl`** (via `makeIndelTSV`): extend output to include all indels with a `zygosity` column (`hom` / `het`) and explicit `ref_allele` / `alt_allele` columns.

**`makeGenomicIndelDb`** (Nextflow process): single table covering both hom and het indels:

```sql
CREATE TABLE genomic_indels (
  strain      TEXT NOT NULL,
  sequence_id TEXT NOT NULL,
  position    INTEGER NOT NULL,
  zygosity    TEXT NOT NULL,   -- 'hom' or 'het'
  ref_allele  TEXT NOT NULL,
  alt_allele  TEXT NOT NULL
);
CREATE INDEX idx_indels ON genomic_indels(strain, sequence_id, position);
```

### CDS extraction — `makeCodingData.jl`

This step produces two SQLite DBs: `codingSequences.db` and `codingIndels.db`.

**`codingSequences.db`** — `coding_sequences(strain, transcript_id, sequence)`

Currently, hom indels are embedded in the consensus FASTA so `extract_cds_sequence` must call `fasta_offset` to shift exon slice coordinates before slicing. Under the new architecture the FASTA is always reference-length, so exon coordinates map directly — `fasta_offset` is removed entirely. Sequences stored here are reference-coordinate slices with SNPs applied and ref bases at all indel positions.

**`codingIndels.db`** — `indels(strain, transcript_id, position, shift_amount)`

`project_indels_to_cds` already translates genomic indel positions to CDS-space coordinates. The schema gains three columns so `processSequenceVariations.jl` can apply the indels directly:

```sql
CREATE TABLE indels (
  strain        TEXT    NOT NULL,
  transcript_id TEXT    NOT NULL,
  position      INTEGER NOT NULL,
  shift_amount  INTEGER NOT NULL,
  zygosity      TEXT    NOT NULL,   -- 'hom' or 'het'
  ref_allele    TEXT    NOT NULL,
  alt_allele    TEXT    NOT NULL
);
```

The coordinate projection logic in `project_indels_to_cds` is otherwise unchanged.

### Product reporting — `processSequenceVariations.jl`

Reads `codingSequences.db` (reference-space CDS slice) and `codingIndels.db` (CDS-space positions + alleles). When a variant record overlaps an indel position:

- **Hom indel:** apply ALT allele to the CDS sequence; translate; emit one product to `product.dat`.
- **Het indel:** produce two sequences — one with REF (the slice as-is) and one with ALT applied; translate each independently; emit both products — same format as het SNP products today.

Strains with no indels in a transcript translate the CDS slice directly, unchanged.

## Affected Files

| File | Change |
|------|--------|
| `bin/makeConsensusFastaFromGvcf.py` | Emit ref bases for all indels (hom and het); `N` reserved for low/no coverage only |
| `bin/findValues.pl` | Output all indels with `zygosity`, `ref_allele`, `alt_allele` columns |
| `modules/mergeExperiments.nf` (`makeGenomicIndelDb`) | Single `genomic_indels` table replacing prior hom-only approach |
| `bin/makeCodingData.jl` | Remove `fasta_offset`; reference-coordinate slicing; add `zygosity`/`ref_allele`/`alt_allele` to `codingIndels.db` |
| `bin/processSequenceVariations.jl` | Load both tracks for het indels; emit both products |

## Out of Scope

- Changes to the SNP or het SNP handling path (unchanged).
- Changes to CNV, windowed coverage, or any workflow other than `processSingleExperiment` → `mergeExperiments`.
- Phasing: multiple het indels in one transcript all receive ALT applied simultaneously — no compound het enumeration.
