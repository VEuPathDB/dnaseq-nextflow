# Per-Sample VCF Publishing + Single Merge Point Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove per-experiment VCF merging from `processSingleExperiment` and publish clean per-sample VCFs instead, so `mergeExperiments` becomes the single merge point across all samples.

**Architecture:** Add a `sanitizeVcf` process to `snp.nf` that strips `%` chars and publishes `{sampleName}.vcf.gz`; wire it into `processSingleExperiment` in place of the old `mergeVcfs` + `makeMergedVariantIndex` pair; update the `mergeExperiments` config glob to match per-sample filenames.

**Tech Stack:** Nextflow DSL2, bcftools, bgzip, tabix (all inside `veupathdb/dnaseqanalysis:1.0.0`)

---

## File Map

| File | Change |
|---|---|
| `modules/snp.nf` | Add `sanitizeVcf`; delete `mergeVcfs` and `makeMergedVariantIndex` |
| `workflows/processSingleExperiment.nf` | Swap imports; rename `combinedVcf` → `filteredVcf`; drop `perSampleVcf`; wire `sanitizeVcf` |
| `nextflow.config` | Update `vcfFiles` glob in `mergeExperiments` profile |

---

## Task 1: Add `sanitizeVcf` to `modules/snp.nf`

**Files:**
- Modify: `modules/snp.nf`

- [ ] **Step 1: Add `sanitizeVcf` process after `makeIndelTSV` (before `mergeVcfs` at line 104)**

Insert the following block at line 103 (after the closing `}` of `makeIndelTSV`):

```groovy
process sanitizeVcf {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir", mode: "copy"

  input:
    tuple val(sampleName), path(vcfGz), path(vcfGzTbi)

  output:
    tuple val(sampleName), path("${sampleName}.vcf.gz"), path("${sampleName}.vcf.gz.tbi")

  script:
    """
    set -euo pipefail
    gunzip -c $vcfGz | sed 's/%//g' | bgzip > ${sampleName}.vcf.gz
    tabix -fp vcf ${sampleName}.vcf.gz
    """

  stub:
    """
    touch ${sampleName}.vcf.gz
    touch ${sampleName}.vcf.gz.tbi
    """
}
```

- [ ] **Step 2: Delete the `mergeVcfs` process (lines 104–130)**

Remove the entire block:

```groovy
process mergeVcfs {
  container 'biocontainers/bcftools:v1.9-1-deb_cv1'

  input:
    val vcfCount
    tuple path(samplevcfzip), path(samplevcfzipindex), val(key)

  output:
    path('result.vcf.gz')

  script:
    """
    set -euo pipefail

    if [ $vcfCount -gt 1 ]; then
        bcftools merge -o result.vcf.gz -O z *.vcf.gz
    else
        cp *.vcf.gz result.vcf.gz
    fi
    """

  stub:
    """
    touch result.vcf.gz
    """

}
```

- [ ] **Step 3: Delete the `makeMergedVariantIndex` process (lines 132–160)**

Remove the entire block:

```groovy
process makeMergedVariantIndex {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir", mode: "copy"

  input:
    path(resultVcfGz)

  output:
    tuple path('result.vcf.gz'), path('result.vcf.gz.tbi')

  script:
    """
    set -euo pipefail
    cp $resultVcfGz hold.vcf.gz
    gunzip hold.vcf.gz
    sed -i 's/\\%//g' hold.vcf
    bgzip hold.vcf
    mv hold.vcf.gz result.vcf.gz
    tabix -fp vcf result.vcf.gz
    """

  stub:
    """
    touch result.vcf.gz
    touch result.vcf.gz.tbi
    """

}
```

- [ ] **Step 4: Verify the file has exactly 6 processes (no orphaned braces)**

```bash
grep "^process " modules/snp.nf
```

Expected:
```
process runFreebayes {
process filterAndSplitVcf {
process makeIndelTSV {
process sanitizeVcf {
process makeCoverageBed {
process makeConsensusFromVcfAndBed {
```

- [ ] **Step 5: Commit**

```bash
git add modules/snp.nf
git commit -m "refactor: add sanitizeVcf process; remove mergeVcfs and makeMergedVariantIndex from snp.nf"
```

---

## Task 2: Update `workflows/processSingleExperiment.nf`

**Files:**
- Modify: `workflows/processSingleExperiment.nf`

Current state of the relevant import lines (lines 23–24):
```groovy
include { mergeVcfs } from '../modules/snp.nf'
include { makeMergedVariantIndex } from '../modules/snp.nf'
```

Current state of the relevant channel/call lines (lines 92–112):
```groovy
    // Extract the per-sample unfiltered VCF (sampleName, vcf.gz, vcf.gz.tbi) for downstream use
    combinedVcf = freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, vcfGz, vcfGzTbi)
    }

    // Feed the indels VCF produced by freebayes directly, bypassing the former filterIndels step
    makeIndelTSV(freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, indelsVcfGz)
    }).collectFile(name: 'indels.tsv', storeDir: params.outputDir)

    // NOTE:  Must ensure the order here is consistent for the vcf files and their indexes;  the lists of paths are each sorted
    mergeVcfsResults = mergeVcfs(combinedVcf.count(), combinedVcf.map{ tuple it[1], it[2], "key" }.groupTuple(by: 2, sort: { a, b -> a <=> b } ))

    makeMergedVariantIndexResults = makeMergedVariantIndex(mergeVcfsResults)

    coverageBedResults = makeCoverageBed(gatkResults.bamTuple)

    perSampleVcf = freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, vcfGz, vcfGzTbi)
    }
    makeConsensusFromVcfAndBed(perSampleVcf.join(coverageBedResults), reorderFastaResults)
```

- [ ] **Step 1: Replace the two old imports with the new one**

Replace:
```groovy
include { mergeVcfs } from '../modules/snp.nf'
include { makeMergedVariantIndex } from '../modules/snp.nf'
```

With:
```groovy
include { sanitizeVcf } from '../modules/snp.nf'
```

- [ ] **Step 2: Rename `combinedVcf` to `filteredVcf` and drop the `perSampleVcf` channel; wire `sanitizeVcf`; remove old merge calls**

Replace the entire block from the `combinedVcf` definition through `makeConsensusFromVcfAndBed(...)`:

```groovy
    // Extract the per-sample unfiltered VCF (sampleName, vcf.gz, vcf.gz.tbi) for downstream use
    combinedVcf = freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, vcfGz, vcfGzTbi)
    }

    // Feed the indels VCF produced by freebayes directly, bypassing the former filterIndels step
    makeIndelTSV(freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, indelsVcfGz)
    }).collectFile(name: 'indels.tsv', storeDir: params.outputDir)

    // NOTE:  Must ensure the order here is consistent for the vcf files and their indexes;  the lists of paths are each sorted
    mergeVcfsResults = mergeVcfs(combinedVcf.count(), combinedVcf.map{ tuple it[1], it[2], "key" }.groupTuple(by: 2, sort: { a, b -> a <=> b } ))

    makeMergedVariantIndexResults = makeMergedVariantIndex(mergeVcfsResults)

    coverageBedResults = makeCoverageBed(gatkResults.bamTuple)

    perSampleVcf = freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, vcfGz, vcfGzTbi)
    }
    makeConsensusFromVcfAndBed(perSampleVcf.join(coverageBedResults), reorderFastaResults)
```

With:
```groovy
    filteredVcf = freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, vcfGz, vcfGzTbi)
    }

    makeIndelTSV(freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, indelsVcfGz)
    }).collectFile(name: 'indels.tsv', storeDir: params.outputDir)

    sanitizeVcfResults = sanitizeVcf(filteredVcf)

    coverageBedResults = makeCoverageBed(gatkResults.bamTuple)

    makeConsensusFromVcfAndBed(sanitizeVcfResults.join(coverageBedResults), reorderFastaResults)
```

- [ ] **Step 3: Verify stub run parses and the process graph looks correct**

```bash
nextflow run main.nf -profile processSingleExperiment -stub-run 2>&1 | head -60
```

Expected: workflow executes stub processes without error, no `Unknown process` or `No such variable` errors.

- [ ] **Step 4: Commit**

```bash
git add workflows/processSingleExperiment.nf
git commit -m "refactor: replace per-experiment VCF merge with per-sample sanitizeVcf in processSingleExperiment"
```

---

## Task 3: Update `nextflow.config` `vcfFiles` glob

**Files:**
- Modify: `nextflow.config`

Current line 55:
```groovy
vcfFiles = "$launchDir/data/merge_setup/**/results/result.vcf.gz"
```

- [ ] **Step 1: Update the glob to match per-sample VCF files**

Replace:
```groovy
vcfFiles                     = "$launchDir/data/merge_setup/**/results/result.vcf.gz"
```

With:
```groovy
vcfFiles                     = "$launchDir/data/merge_setup/**/results/*.vcf.gz"
```

- [ ] **Step 2: Verify stub run for mergeExperiments parses correctly**

The `data/merge_setup` test data already has per-sample VCFs at `data/merge_setup/*/results/*.vcf.gz` (from prior runs). Verify the channel sees them:

```bash
nextflow run main.nf -entry mergeExperiments -profile mergeExperiments -stub-run 2>&1 | head -60
```

Expected: no errors, `mergeVcfs` process shown in the execution graph.

- [ ] **Step 3: Commit**

```bash
git add nextflow.config
git commit -m "config: update mergeExperiments vcfFiles glob to match per-sample VCF filenames"
```

---

## Task 4: Verify end-to-end with test suite

- [ ] **Step 1: Run the test suite**

```bash
nextflow run main.nf -entry runTests -profile tests 2>&1 | tail -30
```

Expected: all Perl `Test2::V0` tests pass, no `FAILED` lines in output.

- [ ] **Step 2: Confirm `result.vcf.gz` is gone from published outputs**

After a stub run of `processSingleExperiment`, verify only per-sample files appear:

```bash
nextflow run main.nf -profile processSingleExperiment -stub-run && ls output/
```

Expected: `*.vcf.gz`, `*.vcf.gz.tbi`, `*.coverage.bed.gz`, `*_consensus.fa.gz`, `indels.tsv` — no `result.vcf.gz`.
