# Design: Replace gVCF with VCF + Coverage BED for Consensus FASTA

**Date:** 2026-05-01  
**Branch:** hand-roll-coverage-per-experiment  
**Scope:** `processSingleExperiment` only — `mergeExperiments` addressed separately

---

## Problem

FreeBayes `--gvcf` produces `<*>` span rows that represent coverage blocks. These cause two classes of bugs:

1. Per-sample: span rows have correctness issues requiring `splitGvcfAtZeroCoverage` as a workaround
2. Cross-sample: `bcftools merge` mishandles span rows in the merged gVCF

The fix is to stop using gVCF entirely for coverage tracking and replace it with an explicit BED file derived directly from the BAM.

---

## Architecture

```
BAM ──────────────────────────────┐
                                   ▼
freebayes (VCF only)          makeCoverageBed
     │                             │
     │ per-sample VCF              │ per-sample coverage.bed
     │                             │
     └──────────┬──────────────────┘
                ▼
    makeConsensusFromVcfAndBed ──► consensus.fa.gz (published)
                                   │
    coverageBedResults.collect()   │
                ▼                  │
    mergeCoverageBeds ─────────────► coverage.tsv (published)
```

---

## Section 1 — `freebayes` changes

- Remove `--gvcf` flag from freebayes invocation
- Remove `${sampleName}.g.vcf.gz` and `${sampleName}.g.vcf.gz.tbi` from output tuple
- Delete `splitGvcfAtZeroCoverage` process definition from `snp.nf`
- Delete `bin/splitGvcfAtZeroCoverage.py`

---

## Section 2 — `makeCoverageBed` (new process in `snp.nf`)

**Input:** `(sampleName, bamFile, bamIndex)`  
**Output:** `(sampleName, "${sampleName}.coverage.bed")`

```bash
bedtools genomecov -ibam $bamFile -bga \
  | awk -v mc=$params.minCoverage '$4 >= mc {print $1"\t"$2"\t"$3}' \
  | bedtools merge \
  > ${sampleName}.coverage.bed
```

Produces a minimal, non-overlapping BED of regions with depth >= `minCoverage`. `bedtools merge` collapses adjacent intervals that both pass the depth threshold.

---

## Section 3 — `makeConsensusFromVcfAndBed` (replaces `makeConsensusFromGvcf`)

**Input:**
- `(sampleName, vcfGz, vcfGzTbi)` — unfiltered per-sample VCF from freebayes
- `(sampleName, coverageBed)` — from `makeCoverageBed`
- `(genomeReorderedFasta, genomeReorderedFastaIndex)`

**New script:** `bin/makeConsensusFastaFromVcfAndBed.py` (replaces `makeConsensusFastaFromGvcf.py`)

### Algorithm

1. Load coverage BED into a per-chrom sorted interval list (bisect for O(log n) lookups)
2. Per chrom, walk VCF records (sparse — variant sites only), tracking `ref_pos = 0`:
   - **Gap** `[ref_pos, variant.POS)`: intersect gap with coverage BED intervals
     - Covered sub-intervals → emit `ref_seq[start:end]`
     - Uncovered sub-intervals → emit `N × length`
   - **Variant site**: check if position falls in a covered BED interval
     - Not covered → emit `N × len(REF)`, advance `ref_pos`
     - Covered → apply variant logic below, advance `ref_pos`
   - After last variant: fill `[ref_pos, chrom_len)` with same gap logic
3. No per-base loop — gap fills use direct string slicing against `ref_seq`

### Variant logic (covered positions)

| Case | Output |
|---|---|
| SNP (all alleles length 1) | IUPAC code for unique alleles |
| Homozygous indel | Actual allele sequence (length may differ from REF) |
| Het indel 0/1 (REF in alleles) | REF sequence — no coordinate shift |
| Het indel 1/2 or both-non-ref complex variant | `X × len(REF)` — same length as REF |
| Uncallable (`.` in GT) | `N × len(REF)` |

**Note on het indels:** `indels.tsv` does not capture het indels. The `X × len(REF)` masking signals ambiguity without introducing a coordinate shift. `len(REF)` is used intentionally — not `len(ALT)`.

---

## Section 4 — `mergeCoverageBeds` (new process, replaces `mergeGvcfs` in `snp.nf`)

**Input:** all per-sample `*.coverage.bed` files collected  
**Output:** `coverage.tsv` published to `outputDir`

```bash
names=$(ls *.coverage.bed | sed 's/\.coverage\.bed//' | tr '\n' ' ')
bedtools multiinter -names $names -i *.coverage.bed > coverage.tsv
```

Output columns: `chrom, start, end, num_samples, sample_list, sample1, sample2, ...`

All regions covered by any sample are reported. Samples with no coverage in a region get `0` in their column. Intervals are split at all sample boundaries.

The old `mergeGvcfs` process definition in `snp.nf` (which produced `coverage.g.vcf.gz`) is removed. The `mergeGvcfs` process in `modules/mergeExperiments.nf` is left untouched (out of scope).

---

## Section 5 — Workflow wiring (`processSingleExperiment.nf`)

**Removed imports:**
- `splitGvcfAtZeroCoverage`
- `makeConsensusFromGvcf`
- `mergeGvcfs` (from `snp.nf`)

**Added imports:**
- `makeCoverageBed`
- `mergeCoverageBeds`
- `makeConsensusFromVcfAndBed`

**Removed calls:**
- `freebayesGvcf` channel extraction
- `splitGvcfResults`
- `makeConsensusFromGvcf(...)` call
- `mergeGvcfs(...)` call

**New wiring:**
```groovy
coverageBedResults = makeCoverageBed(gatkResults.bamTuple)

perSampleVcf = freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
    tuple(sampleName, vcfGz, vcfGzTbi)
}
makeConsensusFromVcfAndBed(perSampleVcf.join(coverageBedResults), reorderFastaResults)

mergeCoverageBeds(coverageBedResults.map { sampleName, bed -> bed }.collect())
```

**Existing maps that must also be updated to 7-element destructuring** (gvcf pair removed):
- `combinedVcf` extraction (~line 92)
- `snpsVcf` extraction (~line 152, inside `if (params.ploidy != 1)` block)

---

## Files Changed

| File | Change |
|---|---|
| `modules/snp.nf` | Remove `freebayes --gvcf`, remove `splitGvcfAtZeroCoverage`, remove `makeConsensusFromGvcf`, remove `mergeGvcfs`; add `makeCoverageBed`, `mergeCoverageBeds`, `makeConsensusFromVcfAndBed` |
| `workflows/processSingleExperiment.nf` | Update imports and wiring as described in Section 5 |
| `bin/makeConsensusFastaFromVcfAndBed.py` | New script |
| `bin/makeConsensusFastaFromGvcf.py` | Delete |
| `bin/splitGvcfAtZeroCoverage.py` | Delete |

---

## Out of Scope

- `mergeExperiments` workflow — still consumes gVCFs; will be addressed in a follow-up
- `bcftoolsMpileupGvcf` process — not used in `processSingleExperiment`, left untouched
