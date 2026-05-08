# Design: Per-Sample VCF Publishing + Single Merge Point

**Date:** 2026-05-06  
**Branch:** hand-roll-coverage-per-experiment  

## Problem

`processSingleExperiment` currently merges all per-sample VCFs into a single `result.vcf.gz` per experiment before publishing. `mergeExperiments` then re-merges those per-experiment VCFs. This is a two-stage merge with redundant intermediate files.

## Goal

- `processSingleExperiment` publishes individual per-sample VCFs (with coverage and consensus)
- All merging happens once, in `mergeExperiments`, across all individual samples from all experiments

## Design

### `modules/snp.nf`

**Add** `sanitizeVcf` process:
- Input: `(sampleName, vcfGz, vcfGzTbi)`
- Strips `%` characters from the VCF (FreeBayes output quirk), bgzips, tabixes
- `publishDir "$params.outputDir"` — publishes `{sampleName}.vcf.gz` and `{sampleName}.vcf.gz.tbi`
- Output: `(sampleName, "{sampleName}.vcf.gz", "{sampleName}.vcf.gz.tbi")`

**Remove** `mergeVcfs` process — no longer called from this module.  
**Remove** `makeMergedVariantIndex` process — replaced by `sanitizeVcf`.

### `workflows/processSingleExperiment.nf`

- Remove imports: `mergeVcfs`, `makeMergedVariantIndex`
- Add import: `sanitizeVcf`
- Rename `combinedVcf` channel to `filteredVcf` (name was a misleading holdover)
- Drop redundant `perSampleVcf` channel (identical mapping to `combinedVcf`)
- Call `sanitizeVcf(filteredVcf)` — publishes clean per-sample VCFs
- Feed `sanitizeVcf` output into `makeConsensusFromVcfAndBed` (replaces raw `perSampleVcf`)
- Remove `mergeVcfsResults` and `makeMergedVariantIndexResults` calls

**Published outputs per sample (unchanged except `result.vcf.gz` → `{sampleName}.vcf.gz`):**
- `{sampleName}.vcf.gz` + `.tbi` (new, from `sanitizeVcf`)
- `{sampleName}.coverage.bed.gz` (unchanged)
- `{sampleName}_consensus.fa.gz` (unchanged)
- `indels.tsv` (unchanged, collectFile)

### `workflows/mergeExperiments.nf`

No logic changes. The existing `allVcfsBranch` single/multiple logic handles both the single-sample edge case and the typical multi-sample path. `mergeVcfs` in `mergeExperiments.nf` re-creates `.tbi` files internally so pre-existing indexes don't need to be passed in.

### `nextflow.config` — `mergeExperiments` profile

Change:
```
vcfFiles = "$launchDir/data/merge_setup/**/results/result.vcf.gz"
```
To:
```
vcfFiles = "$launchDir/data/merge_setup/**/results/*.vcf.gz"
```

## What is NOT changing

- Coverage bed publishing (`makeCoverageBed`) — unchanged
- Consensus publishing (`makeConsensusFromVcfAndBed`) — unchanged (just fed sanitized VCF instead of raw)
- `mergeExperiments` merge logic — unchanged
- `indels.tsv` collection — unchanged
