# Multi-Sample FreeBayes Design

**Date:** 2026-04-30
**Branch:** multi-sample-freebayes
**Status:** Approved for implementation

## Problem

The current pipeline runs FreeBayes independently per sample, then merges per-sample VCFs with `bcftools merge`. This post-hoc merge cannot resolve multi-allelic conflicts across samples correctly and produces merge artifacts that require downstream workarounds.

## Goal

Replace the per-sample FreeBayes → bcftools merge pattern with a joint multi-sample FreeBayes call across all samples in an experiment, using genome chunking for parallelism.

---

## Data Flow

```
GATK BAMs (all samples)
    │
    ├─→ makeRegionBed(fai, chunkSize) → region channel
    │
    │   Per-region (parallelized):
    ├─→ freebayesMultiSample(region, all BAMs, genomeFasta)
    │       → multisample.{region}.g.vcf.gz  (gVCF, all samples)
    │       → multisample.{region}.vcf.gz    (variants only, all samples)
    │
    ├─→ makeMultiSampleZeroCoverageBed(all BAMs, region)
    │       → zero_cov.{region}.bed  (union of all-sample zero-coverage positions)
    │
    ├─→ splitGvcfAtZeroCoverage(multi-sample gVCF chunk, zero_cov.bed)
    │       → split multi-sample gVCF chunk
    │
    │   Gather:
    ├─→ concatMultiSampleVcf(all split gVCF chunks)
    │       → multisample.g.vcf.gz + multisample.vcf.gz (joint, all samples)
    │
    ├─→ makeConsensusFromGvcf(multisample.g.vcf.gz)
    │       → {sampleName}_consensus.fa.gz  (one per sample)
    │
    ├─→ makeIndelTSV(multisample.g.vcf.gz)
    │       → indels.tsv
    │
    └─→ extractSampleVcf(sampleName, multisample.vcf.gz)
            → per-sample VCFs for CNV steps (makeSnpDensity, convertFreebayesToVarscanFormat, etc.)
```

---

## New Processes

### `makeRegionBed`
- **Input:** genome `.fai`, `params.chunkSize`
- **Tool:** `bedtools makewindows -g {fai} -w {chunkSize}`
- **Output:** BED file; emits a channel of region strings (`chrN:start-end`)
- **Config:** New param `params.chunkSize` (e.g. `1000000`)

### `makeMultiSampleZeroCoverageBed`
- **Input:** all BAMs (collected), region string
- **Logic:** Runs `bedtools genomecov -ibam -bga` per sample restricted to the region; takes the union of zero-coverage intervals across all samples
- **Output:** `zero_cov.bed` — shared split boundaries for all samples in this region
- **Key invariant:** A reference block is split at any position where *any* sample has zero coverage

### `concatMultiSampleVcf`
- **Input:** all per-region split multi-sample gVCF chunks (sorted by region)
- **Tools:** `bcftools concat --naive-force` → `bcftools sort` → bgzip + tabix
- **Output:** `multisample.g.vcf.gz` + `multisample.g.vcf.gz.tbi` (gVCF with reference blocks) and `multisample.vcf.gz` + `multisample.vcf.gz.tbi` (variants only, reference blocks stripped)
- **Note:** Any character sanitization (e.g. `%` stripping) done here in one place

### `extractSampleVcf`
- **Input:** joint `multisample.vcf.gz`, `sampleName`
- **Tool:** `bcftools view -s {sampleName} --trim-alt-alleles`
- **Output:** per-sample VCF for CNV downstream steps
- **Runs:** Once per sample

---

## Modified Processes

### `freebayesMultiSample`
- Add `region` input (string, e.g. `chr1:1-1000000`)
- Pass `--region {region}` to freebayes
- Output both gVCF (all records) and variant-only VCF (reference blocks stripped via `bcftools view -e 'ALT[0]="<*>"'`) per region

### `splitGvcfAtZeroCoverage`
- **Remove:** per-sample BAM input and internal `bedtools genomecov` call
- **Add:** `zero_cov.bed` input (pre-computed by `makeMultiSampleZeroCoverageBed`)
- Now operates on multi-sample gVCF chunk, preserving all sample columns
- The `splitGvcfAtZeroCoverage.py` script must be updated to accept `--zero-cov-bed` instead of `--bedgraph`

### `makeConsensusFromGvcf`
- **Change:** takes joint multi-sample gVCF (not a per-sample tuple)
- Iterates over all samples in the gVCF, outputs one `{sampleName}_consensus.fa.gz` per sample
- Single process invocation publishes N consensus files

### `makeIndelTSV`
- **Change:** takes joint multi-sample gVCF instead of per-sample indels VCF
- `findValues.pl` called once per sample from the multi-sample gVCF (loop internally or via `-s` flag if already supported)
- Output: `indels.tsv` (all samples, same format as before)

---

## Removed Processes

| Process | Reason |
|---|---|
| `freebayes` (per-sample) | Replaced by chunked `freebayesMultiSample` |
| `mergeVcfs` | Replaced by `concatMultiSampleVcf` |
| `makeMergedVariantIndex` | Absorbed into `concatMultiSampleVcf` |
| `mergeGvcfs` | Replaced by `concatMultiSampleVcf` on gVCF side |
| `bcftoolsMpileupGvcf` | Unused |

---

## Workflow Changes (`processSingleExperiment.nf`)

### Removed channel wiring
- `freebayesResults` and all downstream maps extracting per-sample VCFs from it
- `mergeVcfs`, `makeMergedVariantIndex`, `mergeGvcfs` calls
- The `freebayesMultiSample` call that currently runs as a parallel extra (lines 124–128)

### New channel wiring
```
regions = makeRegionBed(reorderFastaResults, params.chunkSize)

// Scatter: all BAMs × all regions
regionBams = gatkResults.bamTuple
    .map { sampleName, bam, bai -> tuple(bam, bai) }
    .collect()
    .combine(regions)

multiSampleChunks = freebayesMultiSample(regionBams, reorderFastaResults)

zeroCovBeds = makeMultiSampleZeroCoverageBed(
    gatkResults.bamTuple.map { sampleName, bam, bai -> bam }.collect(),
    regions
)

splitChunks = splitGvcfAtZeroCoverage(
    multiSampleChunks.gvcf.join(zeroCovBeds, by: regionKey),
    reorderFastaResults
)

// Gather
concatResult = concatMultiSampleVcf(splitChunks.collect())

makeConsensusFromGvcf(concatResult.gvcf, reorderFastaResults)
makeIndelTSV(concatResult.gvcf)

// Per-sample extraction for CNV steps
sampleNames = gatkResults.bamTuple.map { sampleName, bam, bai -> sampleName }
perSampleVcf = extractSampleVcf(sampleNames.combine(concatResult.vcf))
```

### CNV steps
`makeSnpDensity`, `convertFreebayesToVarscanFormat`, `getHeterozygousSNPs` receive `perSampleVcf` channel instead of pulling from `freebayesResults`.

---

## Configuration

Add to `nextflow.config` under `processSingleExperiment` profile:

```groovy
chunkSize = 1000000  // 1 Mb windows; tune based on genome size and sample count
```

---

## Key Invariants

1. `splitGvcfAtZeroCoverage` always runs **per-region before concatenation** — never on a whole-genome gVCF
2. Zero-coverage splits are derived from the **union** of all sample coverage bedgraphs — no sample can have a reference block that another sample cannot support
3. The joint multi-sample VCF is the single source of truth for variant calls — no per-sample re-calling or post-hoc merge

---

## Open Questions

- Does `findValues.pl` already handle multi-sample VCF via a loop/flag, or does `makeIndelTSV` need a shell loop over sample names?
- Should `concatMultiSampleVcf` use `--naive-force` (assumes matching headers) or a full `bcftools merge` on the chunks? FreeBayes chunks from the same run should have identical headers, so `--naive-force` should be safe.
