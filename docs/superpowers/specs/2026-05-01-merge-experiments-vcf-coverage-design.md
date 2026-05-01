# Design: Update mergeExperiments to Use VCF + Coverage BED Instead of gVCF

**Date:** 2026-05-01  
**Branch:** hand-roll-coverage-per-experiment  
**Scope:** `mergeExperiments` workflow and `processSequenceVariations.jl`

---

## Problem

`processSingleExperiment` no longer produces gVCF files (removed in the VCF+coverage-BED redesign). It now produces:
- Per-experiment merged VCF (`mergeVcfs` output)
- Per-sample `*.coverage.bed` files (from `makeCoverageBed`, mean depth per interval)
- Per-experiment `consensus.fa.gz` and `indels.tsv` (unchanged)

`mergeExperiments` currently expects per-experiment gVCF files. `processSequenceVariations.jl` uses gVCF REF blocks to determine per-sample coverage at non-variant positions. Both must be updated.

---

## Approach

- Individual per-sample `*.coverage.bed` files are the durable artifact published by `processSingleExperiment`. The `mergeCoverageBeds` process is removed from that workflow — its only purpose was feeding `mergeExperiments`, so it belongs there.
- `mergeExperiments` collects all per-sample bed files across all experiments and merges them into one `coverage.tsv` via `bedtools unionbedg`.
- `processSequenceVariations.jl` reads one `coverage.tsv` with a single file handle, loading one chromosome's intervals at a time to bound memory.
- `bcftools merge` (via existing `mergeVcfs` process) handles VCF merging. `mergeGvcfs` is deleted.
- VCF and `coverage.tsv` must share the same chromosome sort order — guaranteed since both derive from the same reordered reference FASTA.

---

## Section 1 — `processSingleExperiment` changes

**Remove `mergeCoverageBeds`:** Delete the process call from `workflows/processSingleExperiment.nf` and the process definition from `modules/snp.nf`. Remove the corresponding import.

**Publish individual bed files:** Add `publishDir` to `makeCoverageBed` so each `${sampleName}.coverage.bed.gz` is written to `params.outputDir`. The process output is piped through `bgzip` before publishing to save space. `bedtools unionbedg` handles gzipped inputs natively so no decompression step is needed in `mergeExperiments`.

**Dependency on `makeConsensusFastaFromVcfAndBed.py`:** `makeCoverageBed` passes its bed output to `makeConsensusFromVcfAndBed`. Changing the output to `.gz` requires `makeConsensusFastaFromVcfAndBed.py` to detect and open gzipped BED input (e.g., `gzip.open` when filename ends in `.gz`).

No other changes to `processSingleExperiment`.

---

## Section 2 — `processSequenceVariations.jl`

### New argument

`--coverage_file` — single path to the merged `coverage.tsv` produced by `mergeExperiments`.

### Coverage file format

```
chrom    start    end    sample1    sample2    ...
LmjF.01  0        61     30.37      0
LmjF.01  61       14689  30.37      38.90
```

Header line present. Values are mean depth (float); `0` means no coverage.

### Startup

Open the file, read the header to extract sample names and their column indices. Store as `sample_col_index::Dict{String, Int}`. Keep the IO handle open with a peeked-line buffer.

### Per-chromosome coverage loading

The main loop gains a `current_chrom::String` tracker (initially `""`). When the VCF record's chrom differs from `current_chrom`:

1. Set `current_chrom = record.chrom`
2. Call `load_chrom_coverage!(coverage_fh, current_chrom, chrom_coverage)`
   - Advance past any lines for prior chroms
   - Read all lines where `fields[1] == current_chrom`
   - For each sample with `dp > 0`, push `(start, end, dp)` onto the sample's interval vector
   - Stop when the peeked line is a new chrom or EOF
3. Interval vectors are sorted (input is sorted)

`chrom_coverage::Dict{String, Vector{Tuple{Int,Int,Float64}}}` maps sample name → sorted `(start, end, mean_dp)` tuples. Replaced at each chrom transition.

### Coverage lookup

```julia
function get_coverage(chrom_coverage, sample, pos) -> (covered::Bool, mean_dp::Float64)
    intervals = get(chrom_coverage, sample, nothing)
    isnothing(intervals) && return (false, 0.0)
    idx = searchsortedlast(intervals, pos, by = x -> x[1])
    idx == 0 && return (false, 0.0)
    (_, end_, dp) = intervals[idx]
    return pos <= end_ ? (true, dp) : (false, 0.0)
end
```

### Replace `prev_coverage_span`

In `build_variations_from_record` and `handle_variant_record!`:
- Remove `prev_coverage_span::Dict{String, Tuple{String,Int,Int}}` parameter
- Add `chrom_coverage` parameter
- Replace every `prev_coverage_span` lookup with `get_coverage(chrom_coverage, strain, record.pos)`
- DP value comes from `mean_dp` returned by `get_coverage`

### Remove REF block handling

- Delete the `if record.is_ref_block` branch (lines 1641–1661) from the main loop
- Remove `is_ref_block::Bool` and `end_pos::Int` fields from `GVCFRecord`; rename struct to `VCFRecord`
- Remove `peek_end_key` (used only for span-aware REF block sort key)

### Cache drain — no changes

The `cache_key < vcf_start` drain at the top of the main loop already handles positions that fall between variant records, covering the same positions the REF block drain handled.

---

## Section 3 — `modules/mergeExperiments.nf`

**Delete `mergeGvcfs`:** Used `--merge all` for gVCF span rows. No longer needed.

**`mergeVcfs` is already present** in the module (plain `bcftools merge`). It becomes the VCF merge step.

**Add `mergeCoverageBeds` process:**

```groovy
process mergeCoverageBeds {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    path "*.coverage.bed.gz"

  output:
    path 'coverage.tsv'

  script:
    """
    set -euo pipefail
    files=( *.coverage.bed.gz )
    names=( "\${files[@]/.coverage.bed.gz/}" )
    header="chrom\tstart\tend\t\$(IFS='\t'; echo "\${names[*]}")"
    echo -e "\$header" > coverage.tsv
    bedtools unionbedg -names "\${names[@]}" -filler 0 -i "\${files[@]}" >> coverage.tsv
    """

  stub:
    """
    touch coverage.tsv
    """
}
```

**Update `processSeqVars` input and script:**

```groovy
input:
  path vcfFile
  path cacheFile
  path undoneStrainsFile
  val  reference_strain
  path transcriptDb
  path indelDb
  path gtfFile
  path coverageFile    // new

script:
  """
  processSequenceVariations.jl \
    --vcf_file $vcfFile \
    --cache_file $cacheFile \
    --undone_strains_file $undoneStrainsFile \
    --reference_strain $reference_strain \
    --transcript_db $transcriptDb \
    --indel_db $indelDb \
    --gtf_file $gtfFile \
    --coverage_file $coverageFile
  ...
  """
```

---

## Section 4 — `workflows/mergeExperiments.nf`

```groovy
include { mergeCoverageBeds } from '../modules/mergeExperiments.nf'
include { mergeVcfs }         from '../modules/mergeExperiments.nf'   // replaces mergeGvcfs

workflow me {

  take:
    fastas_qch
    vcfs_qch          // renamed from gvcfs_qch
    indels_qch
    coverages_qch     // new: per-sample *.coverage.bed files

  main:

    combinedIndels = indels_qch.collectFile(name: 'indel.tsv')
    genomicIndelDb = makeGenomicIndelDb(combinedIndels)
    allFastas      = fastas_qch.collect()

    allVcfs       = vcfs_qch.collect()
    allVcfsBranch = allVcfs.branch { single: it.size() == 1; multiple: true }
    mergedVcf     = allVcfsBranch.single.map { it[0] }
                      .mix(mergeVcfs(allVcfsBranch.multiple))

    coverageTsv = mergeCoverageBeds(coverages_qch.collect())

    codingData = makeCodingData(allFastas, genomicIndelDb, params.gtfFile, params.genomeFastaFile)

    processSeqVarsResults = processSeqVars(
      mergedVcf,
      params.vcfCacheFile,
      params.undoneStrains,
      params.reference_strain,
      codingData.codingSequencesDb,
      codingData.codingIndelsDb,
      params.gtfFile,
      coverageTsv
    )

    snpEff(processSeqVarsResults.outputVcf, params.gtfFile, params.genomeFastaFile)
}
```

---

## Section 5 — `main.nf` and `nextflow.config`

### `main.nf` (`mergeExperiments` entry point)

```groovy
fastas_qch    = Channel.fromPath(params.relativeConsensusFilePattern)
vcfs_qch      = Channel.fromPath(params.vcfFiles)         // was: gVcfFiles
indels_qch    = Channel.fromPath(params.indelsFiles)
coverages_qch = Channel.fromPath(params.coverageFiles)    // new: *.coverage.bed.gz glob

me(fastas_qch, vcfs_qch, indels_qch, coverages_qch)
```

### `nextflow.config` (`mergeExperiments` profile)

| Old param | New param | Description |
|-----------|-----------|-------------|
| `gVcfFiles` | `vcfFiles` | Glob for per-experiment merged VCFs published by `processSingleExperiment` |
| _(absent)_ | `coverageFiles` | Glob for per-sample `*.coverage.bed.gz` files published by `processSingleExperiment` |

---

## Files Changed

| File | Change |
|------|--------|
| `bin/processSequenceVariations.jl` | Add `--coverage_file`; add per-chrom coverage loading; add `get_coverage` lookup; replace `prev_coverage_span` with `chrom_coverage`; remove REF block branch; rename `GVCFRecord` → `VCFRecord`; remove `peek_end_key` |
| `bin/makeConsensusFastaFromVcfAndBed.py` | Handle gzipped BED input (detect `.gz` extension, use `gzip.open`) |
| `modules/snp.nf` | Remove `mergeCoverageBeds` process; bgzip output in `makeCoverageBed`; add `publishDir` to `makeCoverageBed` |
| `workflows/processSingleExperiment.nf` | Remove `mergeCoverageBeds` import and call |
| `modules/mergeExperiments.nf` | Delete `mergeGvcfs`; add `mergeCoverageBeds`; update `processSeqVars` input/script |
| `workflows/mergeExperiments.nf` | Replace `gvcfs_qch` with `vcfs_qch` + `coverages_qch`; replace `mergeGvcfs` with `mergeVcfs`; add `mergeCoverageBeds` call |
| `main.nf` | Replace `gVcfFiles` with `vcfFiles` + `coverageFiles`; update `me()` call |
| `nextflow.config` | Rename `gVcfFiles` → `vcfFiles`; add `coverageFiles` in `mergeExperiments` profile |

---

## Out of Scope

- `loadSingleExperiment` — no dependency on gVCF format
- Verifying the per-experiment merged VCF is already published by `processSingleExperiment` (add `publishDir` to `mergeVcfs` in `modules/snp.nf` if not)
