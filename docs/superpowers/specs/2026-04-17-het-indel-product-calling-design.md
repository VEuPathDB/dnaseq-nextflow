# Het Indel Product Calling — Design Spec

**Date:** 2026-04-17  
**Status:** Approved

## Problem

Heterozygous indels cannot be represented in a single-sequence consensus FASTA because the two alleles differ in length. The current pipeline masks het indel positions with `N`, the same character used for zero-coverage positions. This means:

1. Downstream CDS extraction at het indel positions yields `N` bases → unknown amino acids.
2. Products cannot be called for het indels, even though het SNPs report both products via IUPAC expansion.
3. No distinction is made between "no coverage" and "covered but heterozygous indel."

## Goal

Report both REF-allele and ALT-allele products for heterozygous indels, consistent with how het SNPs are handled today.

## Design

### Signal layer — `makeConsensusFastaFromGvcf.py`

Change the het indel masking character from `N` to `X`.

- `N` = no coverage (zero depth or below `minCoverage`)
- `X` = position has coverage but carries a heterozygous indel that cannot be collapsed into a single allele

`X` is not a nucleotide IUPAC code and not an amino acid IUPAC code (uppercase `X` is "unknown amino acid" in protein FASTA, but this is a nucleotide FASTA). It is unambiguous in this context.

**Change:** `bin/makeConsensusFastaFromGvcf.py`, line ~106:
```python
# before
segments.append('N' * len(v.REF))
# after
segments.append('X' * len(v.REF))
```

### Data propagation — het indel alleles into the pipeline

`makeCodingData.jl` only has the consensus FASTA and the genomic indels DB. To generate the ALT-track CDS sequence, it needs the REF and ALT allele sequences at each het indel position.

**`findValues.pl`** (via `makeIndelTSV`): extend output to include het indels, with a `is_het` flag and explicit `ref` / `alt` allele columns.

**`makeGenomicIndelDb`** (Nextflow process): add a second table:

```sql
CREATE TABLE genomic_het_indels (
  strain      TEXT NOT NULL,
  sequence_id TEXT NOT NULL,
  position    INTEGER NOT NULL,
  ref_allele  TEXT NOT NULL,
  alt_allele  TEXT NOT NULL
);
CREATE INDEX idx_het ON genomic_het_indels(strain, sequence_id, position);
```

### Dual-track CDS — `makeCodingData.jl`

Extend `coding_sequences` schema with an `allele_track` column:

```sql
CREATE TABLE coding_sequences (
  strain        TEXT NOT NULL,
  transcript_id TEXT NOT NULL,
  allele_track  TEXT NOT NULL DEFAULT 'main',
  sequence      TEXT NOT NULL,
  PRIMARY KEY (strain, transcript_id, allele_track)
);
```

When extracting CDS for a strain:

1. **REF-track (`allele_track = 'main'`):** At `X` positions in the consensus FASTA, substitute the REF allele from the reference genome (no coordinate shift for the het indel). This is the same coordinate space as today.
2. **ALT-track (`allele_track = 'het_alt'`):** Apply the ALT allele from `genomic_het_indels` at those positions. This sequence may differ in length from the REF-track CDS from the het indel position onward.

Only transcripts that overlap at least one het indel position (i.e., contain an `X` in their CDS slice) generate a `het_alt` row. All others store only `allele_track = 'main'`.

### Product reporting — `processSequenceVariations.jl`

When the variant stream encounters a het indel position (identified by `X` in the merged VCF or by the het indel tables):

1. Load both `main` and `het_alt` CDS sequences for the affected `(strain, transcript_id)`.
2. Translate each independently (with existing `expand_codon` / `translate_codon` logic).
3. Emit both products to `product.dat` — same format as het SNP products today.

Strains with no het indels in a transcript continue to use the single `main` sequence path unchanged.

## Affected Files

| File | Change |
|------|--------|
| `bin/makeConsensusFastaFromGvcf.py` | `N` → `X` for het indels |
| `bin/findValues.pl` | Output het indels with ref/alt alleles |
| `modules/mergeExperiments.nf` (`makeGenomicIndelDb`) | Add `genomic_het_indels` table |
| `bin/makeCodingData.jl` | Dual-track CDS generation; schema extension |
| `bin/processSequenceVariations.jl` | Load both tracks; emit both products for het indels |

## Out of Scope

- Changes to the SNP or het SNP handling path (unchanged).
- Changes to CNV, windowed coverage, or any workflow other than `processSingleExperiment` → `mergeExperiments`.
- Phasing: if a transcript has multiple het indels, the `het_alt` track applies all ALT alleles simultaneously (not per-permutation). We do not enumerate compound het combinations.
