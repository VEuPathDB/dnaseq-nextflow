# E2E Test Spec: `mergeExperiments` Workflow

**Date:** 2026-05-08
**Scope:** Spec for automated structural/schema E2E tests of the `mergeExperiments` workflow.
**Input format:** Current format — per-sample `vcf.gz` + `coverage.bed.gz` (lmaj dataset).
**Not covered:** gVCF-based migration branch (spec that branch separately when it lands).

---

## Reference Dataset

| Item | Value |
|---|---|
| Run directory | `/home/jbrestel/dnaseq_test/merge` |
| Strain groups | Mottran (Friedlin_resequence, LV39), WhitePaper (LV39cl5, Seidman751) |
| N strains | 4 |
| Reference strain | `lmajFriedlin` |
| Genome | chr1, haploid |
| GTF | `lmajF_chr1.gtf` |

Input per-strain files (under `input/<group>/results/`):
- `<strain>.vcf.gz` + `<strain>.vcf.gz.tbi` — per-sample FreeBayes VCF
- `<strain>.coverage.bed.gz` — per-base coverage BED
- `<strain>_consensus.fa.gz` — consensus FASTA
- `indels.tsv` — genomic-space indels (one file per group, shared across strains in that group)

---

## Test Approach

Two layers:
1. **Layer 1** — per-output structural invariants for each intermediate and final file
2. **Layer 2** — cross-output consistency checks across all outputs

Tests assert structural correctness only. Exact output values are not pinned (no golden-file comparison). Annotation correctness for `makeCodingData` is covered separately in `docs/qa-makeCodingData-2026-03-18.md`; the checks from that doc are referenced below as promoted-to-automated invariants.

---

## Layer 1: Per-Output Structural Invariants

### `makeGenomicIndelDb` → `genomicIndels.db`

- Valid SQLite file (opens without error)
- Table `genomic_indels` exists with columns: `strain TEXT`, `sequence_id TEXT`, `position INTEGER`, `shift INTEGER`
- Index `idx_genomic_indels` exists on `(sequence_id, strain, position)`
- Row count > 0
- `shift != 0` for all rows (a zero-shift indel is meaningless)

---

### `mergeCoverageBeds` → `coverage.tsv`

- Tab-delimited; header row present as first line
- Header columns: `chrom`, `start`, `end`, followed by one column per input strain
- Total column count == 3 + N_strains
- Strain column names match basenames of input `.coverage.bed.gz` files (strip `.coverage.bed.gz` suffix)
- `start < end` for every data row
- No negative values in any coverage column
- Row count > 0

---

### `mergeVcfs` → `merged.vcf.gz`

- Valid bgzipped VCF (passes `bcftools view` without error)
- Tabix index present (`merged.vcf.gz.tbi`)
- Sample names in VCF header == set of all input strain names
- Every record has at least one non-ref `GT` across all samples
- Record count > 0

Single-sample bypass: when only one input VCF is provided, it is passed through directly (no merge). The same invariants apply.

---

### `makeCodingData` → `codingSequences.db`, `codingIndels.db`

Invariants promoted from `docs/qa-makeCodingData-2026-03-18.md`.

**`codingSequences.db`** — table `coding_sequences(strain, transcript_id, sequence)`:
- Row count == N_strains × N_transcripts_in_GTF (all transcripts present for all strains)
- All strains from input consensus FASTAs are represented
- `sequence` is non-empty for every row
- Every sequence starts with `ATG` (case-insensitive) — Check 1 in QA doc
- Sequence characters match `[ACGTacgtn]+` only
- **Length invariant (Check 3 in QA doc):** for every `(strain, transcript_id)`:
  ```
  len(sequence) == sum(exon_stop - exon_start + 1 for each CDS exon in GTF)
                 + sum(shift for all genomic indels within any exon of this transcript, for this strain)
  ```

**`codingIndels.db`** — table `indels(strain, transcript_id, position, shift_amount)`:
- All strains represented
- `shift_amount != 0` for all rows
- `position >= 1` for all rows (1-based CDS coordinate)

---

### `processSeqVars` → `output.vcf.gz`, `variationFeature.dat`, `allele.dat`, `product.dat`

**`output.vcf.gz`:**
- Valid bgzipped VCF; tabix index present (`output.vcf.gz.tbi`)
- VCF header contains `##INFO=<ID=CANN`
- VCF header contains `##FORMAT=<ID=CA`
- VCF header contains `##FORMAT=<ID=DFS`
- Every record has a non-empty `CANN` INFO field (presence required; ref-only values are valid — see CANN note below)
- Sample names match those in `merged.vcf.gz`
- Record count > 0

**CANN annotation note:** A ref-only or absent `CA` per-sample value is valid when: `DFS=1` (the sample carries an upstream frameshift at this position), the sample has no coverage at the site, or the site falls within a het indel region. No automated assertion is made on CANN content beyond field presence.

**`variationFeature.dat`** — 18 tab-delimited columns, no header:

| Col | Name | Invariant |
|-----|------|-----------|
| 1 | location | positive integer |
| 2 | transcript_id | string (may be empty for non-coding) |
| 3 | seq_id | non-empty string |
| 4 | reference_strain | equals configured `reference_strain` param for every row |
| 5 | ref_allele | non-empty string |
| 6 | has_nonsynonymous_allele | 0 or 1 |
| 7 | major_allele | non-empty string |
| 8 | minor_allele | string (may be empty) |
| 9 | major_allele_count | positive integer |
| 10 | minor_allele_count | integer or empty |
| 11 | major_product | string (may be empty) |
| 12 | minor_product | string (may be empty) |
| 13 | distinct_strain_count | 1 ≤ value ≤ N_strains |
| 14 | distinct_allele_count | ≥ 1 |
| 15 | is_coding | 0 or 1 |
| 16 | total_allele_count | ≥ distinct_strain_count |
| 17 | has_stop_codon | 0 or 1 |
| 18 | ref_codon | string (may be empty) |

Additional row-level invariants:
- `is_coding=1` rows must have non-empty `transcript_id`
- `has_nonsynonymous_allele=1` implies `is_coding=1`
- Row count > 0

**`allele.dat`** — 5 tab-delimited columns, no header; only emitted for coding variants (`is_coding=1`):

| Col | Name | Invariant |
|-----|------|-----------|
| 1 | allele | non-empty string |
| 2 | distinct_strain_count | positive integer |
| 3 | allele_count | ≥ distinct_strain_count |
| 4 | avg_coverage | float ≥ 0.0, 2 decimal places |
| 5 | avg_percent | float 0.0–100.0, 2 decimal places |

**`product.dat`** — 8 tab-delimited columns, no header; only emitted for coding variants with non-empty codons:

| Col | Name | Invariant |
|-----|------|-----------|
| 1 | expanded_codon | 3-character string matching `[ACGT]{3}` |
| 2 | pos_in_codon_val | 1, 2, or 3 |
| 3 | transcript_id | non-empty string |
| 4 | product_count | positive integer |
| 5 | amino_acid | single character or `*` (stop codon) |
| 6 | pos_in_cds | positive integer |
| 7 | pos_in_codon_val (dup) | 1, 2, or 3 |
| 8 | downstream_of_frameshift | 0 or 1 |

---

### `snpEff` → `merged.ann.vcf.gz`

- Valid bgzipped VCF; tabix index present (`merged.ann.vcf.gz.tbi`)
- VCF header contains `##INFO=<ID=ANN` (snpEff annotation)
- VCF header contains `##INFO=<ID=CANN` (carried through from `output.vcf.gz`)
- VCF header contains `##FORMAT=<ID=CA` and `##FORMAT=<ID=DFS` (carried through)
- Every record has a non-empty `ANN` INFO field
- Every record has a non-empty `CANN` INFO field
- Record count == record count in `output.vcf.gz` (snpEff annotates; it does not filter records)

---

## Layer 2: Cross-Output Consistency

### 1. Strain name consistency

- VCF sample names in `merged.vcf.gz` header == set of strain names derived from input file basenames (strip `.vcf.gz` suffix)
- Coverage TSV strain column names (cols 4+) == VCF sample names in `merged.vcf.gz` (same set; order may differ)
- Distinct `strain` values in `genomic_indels` table == distinct `strain` values in `coding_sequences` table == VCF sample names

### 2. Variant count pipeline integrity

- Record count in `output.vcf.gz` == record count in `merged.ann.vcf.gz`
- `variationFeature.dat` row count > 0 and ≤ `output.vcf.gz` record count × max_transcripts_per_position

  The relationship is not 1-to-1: a VCF record that overlaps N transcripts produces N rows; a cache-passthrough record where all current strains share the same allele (`has_variation` = false) produces 0 rows. No exact equality is asserted.

### 3. Coding variant cross-file consistency

- Every row in `allele.dat` corresponds to a `variationFeature.dat` row where `is_coding=1`
- Every `transcript_id` in `product.dat` also appears in `codingSequences.db`

### 4. Reference strain integrity

- The configured `reference_strain` appears as a sample name in `merged.vcf.gz`
- Column 4 of `variationFeature.dat` equals `reference_strain` for every row

### 5. CANN annotation propagation

- Every record in `output.vcf.gz` has a `CANN` INFO field present
- Every record in `merged.ann.vcf.gz` has both `CANN` and `ANN` INFO fields present
- Content of CANN is not asserted beyond presence (see CANN annotation note in processSeqVars section)

---

## Out of Scope

- gVCF input format (hand-roll-coverage-per-experiment branch) — spec separately
- Golden-file / exact-value comparison of any output
- `makeCodingData` annotation correctness beyond structural invariants — covered by `docs/qa-makeCodingData-2026-03-18.md`
- `loadSingleExperiment` (database loading) — separate workflow
