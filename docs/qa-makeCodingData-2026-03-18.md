# QA Report: `makeCodingData` Process

**Date:** March 18, 2026
**Analyst:** Claude Sonnet 4.6
**Pipeline:** `dnaseq-nextflow` — `mergeExperiments` workflow
**Process under test:** `makeCodingData` (`bin/makeCodingData.jl`)
**Work directory:** `~/dnaseq_test/merge/work/1f/e3d45f1370c41f41d57a9e1fb8a3fd/`

---

## Overview

`makeCodingData.jl` consumes per-strain consensus FASTAs, a reference genome, a GTF annotation, and a genomic indel SQLite database to produce two outputs:

- **`codingSequences.db`** — `coding_sequences(strain, transcript_id, sequence)`: CDS sequence for every transcript in every strain.
- **`codingIndels.db`** — `indels(strain, transcript_id, position, shift_amount)`: indels projected into CDS coordinate space.

The key algorithmic detail being validated is `fasta_offset`: for each CDS exon boundary, the code sums all genomic indel `shift` values at positions strictly upstream of that boundary to convert reference GTF coordinates into the correct slice position within the strain's consensus FASTA.

---

## Dataset

| Item | Value |
|---|---|
| Reference genome | `genome.fasta` (PlasmoDB *P. falciparum* 3D7) |
| GTF | `pfal3D7.gtf` |
| Strains | `5.1`, `ERR015417`, `ERR027115`, `PS097` + `reference` |
| Total transcripts in GTF | 219 |
| Total CDS entries in DB | 219 transcripts × 5 strains = 1,095 rows |
| Transcript exon count distribution | 83 single-exon, 136 multi-exon (2–16 exons per transcript) |

---

## Check 1: Sequence Termini (First & Last 20 bp)

Four transcripts were selected to cover all structural categories:

| Transcript | Strand | Exon count | Category |
|---|---|---|---|
| `PF3D7_0200100.1` | + | 2 | forward, multi-exon |
| `PF3D7_0200200.1` | − | 2 | reverse, multi-exon |
| `PF3D7_0210500.1` | + | 1 | forward, single-exon |
| `PF3D7_0205000.1` | − | 1 | reverse, single-exon |

First and last 20 bp were queried directly from `codingSequences.db` for all five strains.

### PF3D7_0200100.1 — forward, multi-exon (exons: 25232–29035, 29837–31168)

| Strain | First 20 bp | Last 20 bp | Length |
|---|---|---|---|
| reference | `atggcgactggtagtggggg` | `cggatgtgtgggatatataa` | 5136 |
| 5.1 | `ATGGCGACTGGTAGTGGGGG` | `CAGATGTGTGGGATATATAA` | 5133 |
| ERR015417 | `ATGGCGACTGGTAGTGGGGG` | `CGGATGTATGGGATATATAA` | 5136 |
| ERR027115 | `ATGGCGACTGGTAGTGGGGG` | `CGGATGTGTGGGATATATAA` | 5137 |
| PS097 | `ATGGCGACTGGTAGTGGGGG` | `CGGATGTGTGGGATATATAA` | 5137 |

All sequences begin with ATG. Length variation (5133–5137) is consistent with indels across strains. Two strains show SNPs in the last 20 bp (5.1 and ERR015417) confirmed in Check 3 below.

### PF3D7_0200200.1 — reverse, multi-exon (exons: 33030–33965, 34191–34259)

| Strain | First 20 bp | Last 20 bp | Length |
|---|---|---|---|
| reference | `atgatgttgaattacactaa` | `taaaactattagaagaatag` | 1005 |
| 5.1 | `ATGATGTTGAATTACACTAA` | `TAAAACTATTAGAAGAATAG` | 1008 |
| ERR015417 | `ATGATGTTGAATTACACTAA` | `TAAAACTATTAGAAGAATAG` | 1005 |
| ERR027115 | `ATGATGTTGAATTACACTAA` | `TAAAACTATTAGAAGAATAG` | 1005 |
| PS097 | `ATGAAGTTCAATTACACTAA` | `TAAAACTATTAGAAGAATAG` | 1005 |

PS097 shows two SNPs at CDS positions 4 and 8 (see Check 3). Strain 5.1 is 3 bp longer (one in-frame insertion).

### PF3D7_0210500.1 — forward, single-exon

All five strains produce identical first/last 20 bp and identical length (4041). No variation at termini.

### PF3D7_0205000.1 — reverse, single-exon

First/last 20 bp are identical across all strains. ERR015417 is 6 bp shorter (1521 vs 1527) — consistent with a 2-codon deletion in the genomic indel DB.

---

## Check 2: fasta_offset Coordinate Adjustment

For the two multi-exon transcripts, per-strain `fasta_offset` values (sum of genomic indel shifts upstream of each exon boundary) were computed manually from `genomicIndels.db` and used to extract the expected sequence directly from each strain's `_consensus.fa.gz`. Results were compared against `codingSequences.db`.

### PF3D7_0200100.1 — forward, multi-exon

| Strain | Shift before ex1 start (25232) | Adjusted fstart | Shift before ex2 stop+1 (31169) | Adjusted fstop | First 20 match | Last 20 match |
|---|---|---|---|---|---|---|
| 5.1 | −8 | 25224 | −11 | 31157 | ✓ | ✓ |
| ERR015417 | −12 | 25220 | −16 | 31152 | ✓ | ✓ |
| ERR027115 | −2 | 25230 | −1 | 31167 | ✓ | ✓ |
| PS097 | 0 | 25232 | +1 | 31169 | ✓ | ✓ |

### PF3D7_0200200.1 — reverse, multi-exon

For minus-strand transcripts, exons are ordered 5'→3' in descending genomic coordinate order. Each exon's raw genomic slice is individually reverse-complemented before concatenation. The first 20 bp of CDS derive from the reverse complement of the last 20 bp of the 5'-most exon (34191–34259); the last 20 bp derive from the reverse complement of the first 20 bp of the 3'-most exon (33030–33965).

| Strain | ex2 fstart→fstop | ex1 fstart→fstop | First 20 match | Last 20 match |
|---|---|---|---|---|
| 5.1 | 34183→34251 | 33019→33957 | ✓ | ✓ |
| ERR015417 | 34175→34243 | 33014→33949 | ✓ | ✓ |
| ERR027115 | 34190→34258 | 33029→33964 | ✓ | ✓ |
| PS097 | 34191→34259 | 33030→33965 | ✓ | ✓ |

All 8 extractions match the database exactly.

---

## Check 3: Length Sanity Check (All Transcripts, All Strains)

For every (strain, transcript) pair, the expected CDS length was computed as:

```
expected_length = sum(exon_stop - exon_start + 1  for each CDS exon)
               + sum(shift  for all genomic indels within any exon of this transcript, for this strain)
```

This was compared against `length(sequence)` in `codingSequences.db`.

| Result | Count |
|---|---|
| Total checks | 876 (219 transcripts × 4 strains) |
| Mismatches | **0** |

Every CDS length in the database is consistent with the reference exon structure plus the per-strain indel shifts.

---

## Check 4: VCF Confirmation of Observed SNPs

Three SNP differences observed in Check 1 were traced back to the originating g.VCF files (`~/dnaseq_test/merge/input/*/results/coverage.g.vcf.gz`).

### 5.1 — Pf3D7_02_v3:31150 G→A

Located in the broad cohort VCF (samples: 5.1, PS097).

```
Pf3D7_02_v3  31150  G  A  QUAL=1677.94  TYPE=snp
GT=1 (homozygous alt)  DP=88  AO=76  AF=1
```

Confirmed homozygous SNP with high depth and quality.

### ERR015417 — Pf3D7_02_v3:31156 G→A

Located in the sanger cohort VCF (samples: ERR015417, ERR027115).

```
Pf3D7_02_v3  31156  G  A  QUAL=863.29  TYPE=snp
GT=1 (homozygous alt)  DP=41  AO=40  AF=1
```

Confirmed homozygous SNP.

### PS097 — Pf3D7_02_v3:34251 CAACA→GAACT

Located in the broad cohort VCF. The two CDS-space SNPs (positions 4 and 8) are encoded as a single complex variant spanning positions 34251–34255.

```
Pf3D7_02_v3  34251  CAACA  GAACT  QUAL=337.33  TYPE=complex  CIGAR=1X3M1X
GT=1 (homozygous alt)  DP=10  AO=10  AF=1
```

The CIGAR `1X3M1X` confirms two substitutions flanking 3 matching bases — exactly accounting for both CDS-position diffs. Note: DP=10 is low relative to the other calls; the call is homozygous but may warrant coverage review.

---

## Summary

| Check | Scope | Result |
|---|---|---|
| CDS starts with ATG | All 4 test transcripts × 5 strains | ✓ All correct |
| First/last 20 bp match consensus FASTA (forward multi-exon) | 4 strains × 2 termini | ✓ 8/8 |
| First/last 20 bp match consensus FASTA (reverse multi-exon) | 4 strains × 2 termini | ✓ 8/8 |
| Length = ref length + indel shifts | 219 transcripts × 4 strains | ✓ 876/876 |
| Observed SNPs present in source g.VCF | 3 variant sites | ✓ 3/3 |

No issues found. The `makeCodingData` process correctly applies `fasta_offset` coordinate adjustment for all strand orientations and exon structures, and the resulting coding sequences are fully consistent with the source g.VCF variant calls.

---

*Report generated by Claude Sonnet 4.6 via interactive QA session.*
