# mergeExperiments VCF + Coverage BED Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace gVCF inputs in `mergeExperiments` with regular per-sample VCFs + per-sample `*.coverage.bed.gz` files, and update `processSequenceVariations.jl` to stream coverage from a merged `coverage.tsv` instead of gVCF REF blocks.

**Architecture:** `processSingleExperiment` publishes bgzipped per-sample bed files; `mergeExperiments` merges all bed files into one `coverage.tsv` via `bedtools unionbedg`, feeds that single file to `processSequenceVariations.jl`, which loads one chromosome's coverage intervals at a time for bounded memory use.

**Tech Stack:** Nextflow DSL2, Julia 1.10 (via `veupathdb/dnaseqanalysis:1.0.0`), bedtools, bcftools, Python 3

**Spec:** `docs/superpowers/specs/2026-05-01-merge-experiments-vcf-coverage-design.md`

---

### Task 1: bgzip `makeCoverageBed` output and publish

**Files:**
- Modify: `modules/snp.nf` — `makeCoverageBed` process

- [ ] **Step 1: Update `makeCoverageBed` to bgzip output and add `publishDir`**

Replace the current `makeCoverageBed` process in `modules/snp.nf`:

```groovy
process makeCoverageBed {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir", mode: "copy"

  input:
    tuple val(sampleName), path(bamFile), path(bamIndex)

  output:
    tuple val(sampleName), path("${sampleName}.coverage.bed.gz")

  script:
    """
    set -euo pipefail
    bedtools genomecov -ibam $bamFile -bga \\
      | awk -v mc=$params.minCoverage '\$4 >= mc' \\
      | bedtools merge -c 4 -o mean \\
      | bgzip > ${sampleName}.coverage.bed.gz
    """

  stub:
    """
    touch ${sampleName}.coverage.bed.gz
    """
}
```

- [ ] **Step 2: Verify stub run passes**

```bash
nextflow run main.nf -profile processSingleExperiment -stub
```

Expected: pipeline completes without errors. `makeCoverageBed` stage present; no reference to the old uncompressed bed filename.

- [ ] **Step 3: Commit**

```bash
git add modules/snp.nf
git commit -m "feat: bgzip makeCoverageBed output and publish to outputDir"
```

---

### Task 2: Update `makeConsensusFastaFromVcfAndBed.py` for gzipped BED input

**Files:**
- Modify: `bin/makeConsensusFastaFromVcfAndBed.py` — `load_coverage_bed` function

- [ ] **Step 1: Replace `load_coverage_bed` with a gzip-aware version**

Replace the existing `load_coverage_bed` function body:

```python
def load_coverage_bed(bed_path):
    """
    Parse a BED file (plain or gzipped) and return a dict mapping
    chrom -> list of (start, end) intervals (0-based, half-open).
    """
    import gzip as _gzip
    opener = _gzip.open if bed_path.endswith('.gz') else open
    result = defaultdict(list)
    with opener(bed_path, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            result[chrom].append((start, end))
    return dict(result)
```

- [ ] **Step 2: Verify no syntax errors**

```bash
python3 -c "import ast; ast.parse(open('bin/makeConsensusFastaFromVcfAndBed.py').read()); print('OK')"
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add bin/makeConsensusFastaFromVcfAndBed.py
git commit -m "fix: handle gzipped BED input in makeConsensusFastaFromVcfAndBed.py"
```

---

### Task 3: Remove `mergeCoverageBeds` from `processSingleExperiment`

**Files:**
- Modify: `modules/snp.nf` — delete `mergeCoverageBeds` process
- Modify: `workflows/processSingleExperiment.nf` — remove import and call

- [ ] **Step 1: Delete the `mergeCoverageBeds` process from `modules/snp.nf`**

Delete the entire `process mergeCoverageBeds { ... }` block (the one with `val sampleNames` / `path bedFiles` inputs and `bedtools unionbedg` script).

- [ ] **Step 2: Remove `mergeCoverageBeds` import from `workflows/processSingleExperiment.nf`**

Delete:
```groovy
include { mergeCoverageBeds } from '../modules/snp.nf'
```

- [ ] **Step 3: Remove `mergeCoverageBeds` call from `workflows/processSingleExperiment.nf`**

Delete:
```groovy
    mergeCoverageBeds(
        coverageBedResults.map { sampleName, bed -> sampleName }.collect(),
        coverageBedResults.map { sampleName, bed -> bed }.collect()
    )
```

- [ ] **Step 4: Verify stub run passes**

```bash
nextflow run main.nf -profile processSingleExperiment -stub
```

Expected: completes without errors. No `mergeCoverageBeds` stage in the output.

- [ ] **Step 5: Commit**

```bash
git add modules/snp.nf workflows/processSingleExperiment.nf
git commit -m "refactor: remove mergeCoverageBeds from processSingleExperiment; bed files published individually"
```

---

### Task 4: Update `modules/mergeExperiments.nf` — delete `mergeGvcfs`, add `mergeCoverageBeds`, update `processSeqVars`

**Files:**
- Modify: `modules/mergeExperiments.nf`

- [ ] **Step 1: Delete the `mergeGvcfs` process**

Delete the entire `process mergeGvcfs { ... }` block (uses `bcftools merge --merge all`).

- [ ] **Step 2: Add `mergeCoverageBeds` process after `mergeVcfs`**

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

- [ ] **Step 3: Add `path coverageFile` as the last input of `processSeqVars`**

The updated input block for `processSeqVars`:

```groovy
  input:
    path vcfFile
    path cacheFile
    path undoneStrainsFile
    val  reference_strain
    path transcriptDb
    path indelDb
    path gtfFile
    path coverageFile
```

- [ ] **Step 4: Add `--coverage_file $coverageFile` to the `processSeqVars` script**

The updated `processSequenceVariations.jl` invocation:

```bash
    processSequenceVariations.jl \\
      --vcf_file $vcfFile \\
      --cache_file $cacheFile \\
      --undone_strains_file $undoneStrainsFile \\
      --reference_strain $reference_strain \\
      --transcript_db $transcriptDb \\
      --indel_db $indelDb \\
      --gtf_file $gtfFile \\
      --coverage_file $coverageFile
```

- [ ] **Step 5: Commit**

```bash
git add modules/mergeExperiments.nf
git commit -m "feat: replace mergeGvcfs with mergeCoverageBeds; wire coverage_file into processSeqVars"
```

---

### Task 5: Update `workflows/mergeExperiments.nf`

**Files:**
- Modify: `workflows/mergeExperiments.nf`

- [ ] **Step 1: Replace the entire file**

```groovy
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { mergeVcfs }          from '../modules/mergeExperiments.nf'
include { mergeCoverageBeds }  from '../modules/mergeExperiments.nf'
include { makeGenomicIndelDb } from '../modules/mergeExperiments.nf'
include { makeCodingData }     from '../modules/mergeExperiments.nf'
include { processSeqVars }     from '../modules/mergeExperiments.nf'
include { snpEff }             from '../modules/mergeExperiments.nf'

workflow me {

  take:
    fastas_qch
    vcfs_qch
    indels_qch
    coverages_qch

  main:

    combinedIndels = indels_qch.collectFile(name: 'indel.tsv')

    genomicIndelDb = makeGenomicIndelDb(combinedIndels)

    allFastas = fastas_qch.collect()

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

- [ ] **Step 2: Commit**

```bash
git add workflows/mergeExperiments.nf
git commit -m "feat: update mergeExperiments workflow — vcfs+coverages_qch, mergeCoverageBeds, mergeVcfs"
```

---

### Task 6: Update `main.nf` and `nextflow.config`

**Files:**
- Modify: `main.nf` — `mergeExperiments` entry point
- Modify: `nextflow.config` — `mergeExperiments` profile params

- [ ] **Step 1: Replace the `workflow mergeExperiments` block in `main.nf`**

Current (lines 81-87):
```groovy
workflow mergeExperiments {
  fastas_qch  = Channel.fromPath(params.relativeConsensusFilePattern)
  gvcfs_qch   = Channel.fromPath(params.gVcfFiles)
  indels_qch  = Channel.fromPath(params.indelsFiles)

  me(fastas_qch, gvcfs_qch, indels_qch)
}
```

Replace with:
```groovy
workflow mergeExperiments {
  fastas_qch    = Channel.fromPath(params.relativeConsensusFilePattern)
  vcfs_qch      = Channel.fromPath(params.vcfFiles)
  indels_qch    = Channel.fromPath(params.indelsFiles)
  coverages_qch = Channel.fromPath(params.coverageFiles)

  me(fastas_qch, vcfs_qch, indels_qch, coverages_qch)
}
```

- [ ] **Step 2: Update `nextflow.config` `mergeExperiments` profile**

Replace:
```groovy
      gVcfFiles = "$launchDir/data/merge_setup/**/results/coverage.g.vcf.gz"
```
with:
```groovy
      vcfFiles      = "$launchDir/data/merge_setup/**/results/result.vcf.gz"
      coverageFiles = "$launchDir/data/merge_setup/**/results/*.coverage.bed.gz"
```

Note: `mergeVcfs` in `modules/snp.nf` publishes its output as `result.vcf.gz`.

- [ ] **Step 3: Verify stub run of `mergeExperiments` passes**

```bash
nextflow run main.nf -entry mergeExperiments -profile mergeExperiments -stub
```

Expected: completes without errors. Stages: `mergeVcfs`, `mergeCoverageBeds`, `makeGenomicIndelDb`, `makeCodingData`, `processSeqVars`, `snpEff`. `mergeGvcfs` does NOT appear.

- [ ] **Step 4: Commit**

```bash
git add main.nf nextflow.config
git commit -m "feat: update main.nf and nextflow.config for vcfFiles/coverageFiles params"
```

---

### Task 7: Add coverage data structures to `processSequenceVariations.jl`

**Files:**
- Modify: `bin/processSequenceVariations.jl` — insert new `Coverage I/O` section after `parse_format_field` (~line 562)

- [ ] **Step 1: Insert the new `Coverage I/O` section**

After the `parse_format_field` function, insert:

```julia
# ---------------------------------------------------------------------------
# Coverage I/O (coverage.tsv produced by mergeCoverageBeds)
# ---------------------------------------------------------------------------

mutable struct CoverageFileHandle
    fh::IO
    sample_cols::Vector{Tuple{Int, String}}   # (col_index_1based, sample_name)
    peeked::Union{String, Nothing}
    exhausted::Bool
end

"""
    open_coverage_file(path) -> CoverageFileHandle

Opens coverage.tsv, reads the header to build a sample→column-index mapping,
and buffers the first data line.
"""
function open_coverage_file(path::String)::CoverageFileHandle
    fh = open(path, "r")
    header = readline(fh)
    fields = split(header, '\t')
    # Columns 4+ are sample names (1-based indexing: column 4 = index 4)
    sample_cols = Tuple{Int, String}[(i, String(fields[i])) for i in 4:length(fields)]
    first_line  = eof(fh) ? nothing : readline(fh)
    CoverageFileHandle(fh, sample_cols, first_line, first_line === nothing)
end

"""
    load_chrom_coverage!(cfh, chrom, chrom_rank, chrom_coverage)

Advances cfh past any lines for chromosomes that sort before `chrom` (by chrom_rank),
then reads all intervals for `chrom` into `chrom_coverage`, replacing prior contents.
Interval vectors are already sorted since coverage.tsv is position-sorted.
"""
function load_chrom_coverage!(
    cfh::CoverageFileHandle,
    chrom::String,
    chrom_rank::Dict{String, Int},
    chrom_coverage::Dict{String, Vector{Tuple{Int, Int, Float64}}}
)
    empty!(chrom_coverage)
    cfh.exhausted && return

    target_rank = get(chrom_rank, chrom, typemax(Int))

    # Advance past chromosomes that sort before the target
    while !cfh.exhausted
        fields    = split(cfh.peeked, '\t')
        line_rank = get(chrom_rank, String(fields[1]), typemax(Int))
        line_rank >= target_rank && break
        cfh.peeked    = eof(cfh.fh) ? nothing : readline(cfh.fh)
        cfh.exhausted = cfh.peeked === nothing
    end

    # Read all lines for this chrom
    while !cfh.exhausted
        fields = split(cfh.peeked, '\t')
        String(fields[1]) != chrom && break

        start_pos = parse(Int, fields[2])
        end_pos   = parse(Int, fields[3])

        for (col_idx, sample) in cfh.sample_cols
            col_idx > length(fields) && continue
            dp = parse(Float64, String(fields[col_idx]))
            if dp > 0.0
                push!(get!(chrom_coverage, sample, Tuple{Int,Int,Float64}[]),
                      (start_pos, end_pos, dp))
            end
        end

        cfh.peeked    = eof(cfh.fh) ? nothing : readline(cfh.fh)
        cfh.exhausted = cfh.peeked === nothing
    end
end

"""
    get_coverage(chrom_coverage, sample, pos) -> (covered, mean_dp)

Binary search for coverage at 0-based `pos`. Returns (false, 0.0) if not covered.
coverage.tsv uses 0-based half-open intervals [start, end) matching BED convention.
Pass VCF positions as `record.pos - 1` to convert from 1-based to 0-based.
"""
function get_coverage(
    chrom_coverage::Dict{String, Vector{Tuple{Int, Int, Float64}}},
    sample::String,
    pos::Int
)::Tuple{Bool, Float64}
    intervals = get(chrom_coverage, sample, nothing)
    (isnothing(intervals) || isempty(intervals)) && return (false, 0.0)
    idx = searchsortedlast(intervals, pos, by = x -> x[1])
    idx == 0 && return (false, 0.0)
    (_, end_, dp) = intervals[idx]
    return pos < end_ ? (true, dp) : (false, 0.0)
end
```

- [ ] **Step 2: Verify syntax using the Docker image**

```bash
docker run --rm -v $(pwd):/work -w /work veupathdb/dnaseqanalysis:1.0.0 \
  julia --compile=min -e 'include("bin/processSequenceVariations.jl")' 2>&1 | head -20
```

Expected: no parse or compilation errors (a runtime error about missing `--vcf_file` arg is OK).

- [ ] **Step 3: Commit**

```bash
git add bin/processSequenceVariations.jl
git commit -m "feat: add CoverageFileHandle, load_chrom_coverage!, get_coverage to processSequenceVariations.jl"
```

---

### Task 8: Remove REF block infrastructure — rename `GVCFRecord`, remove `peek_end_key`

**Files:**
- Modify: `bin/processSequenceVariations.jl`

- [ ] **Step 1: Replace `GVCFRecord` struct with `VCFRecord` (lines 358-368)**

```julia
struct VCFRecord
    chrom::String
    pos::Int
    ref::String
    alts::Vector{String}
    info::String
    format_keys::Vector{String}
    sample_data::Vector{String}   # raw per-sample FORMAT strings
end
```

- [ ] **Step 2: Replace `parse_gvcf_record` with `parse_vcf_record` (lines 538-559)**

```julia
"""
    parse_vcf_record(line, n_samples) -> VCFRecord
"""
function parse_vcf_record(line::String, n_samples::Int)::VCFRecord
    fields = split(line, '\t')
    chrom = String(fields[1])
    pos   = parse(Int, fields[2])
    ref   = String(fields[4])
    alts  = String[String(a) for a in split(fields[5], ',')]
    info  = String(fields[8])
    fmt   = String(fields[9])
    format_keys = String[String(k) for k in split(fmt, ':')]
    sample_data = String[String(fields[9+i]) for i in 1:n_samples if 9+i <= length(fields)]
    VCFRecord(chrom, pos, ref, alts, info, format_keys, sample_data)
end
```

- [ ] **Step 3: Rename `parse_gvcf_header` → `parse_vcf_header` and `open_gvcf_peeked` → `open_vcf_peeked`**

Replace both functions (lines 493-533):

```julia
"""
    parse_vcf_header(io) -> (all_strains, chrom_rank, info_headers)

Reads ## meta lines, builds chrom_rank from ##contig lines, extracts sample
names from #CHROM line. Leaves io positioned at first data line.
"""
function parse_vcf_header(io::IO)
    chrom_rank   = Dict{String,Int}()
    all_strains  = String[]
    info_headers = String[]
    contig_count = 0

    for line in eachline(io)
        if startswith(line, "##INFO")
            push!(info_headers, line)
        elseif startswith(line, "##contig")
            m = match(r"##contig=<ID=([^,>]+)", line)
            if !isnothing(m)
                contig_count += 1
                chrom_rank[m.captures[1]] = contig_count
            end
        elseif startswith(line, "#CHROM")
            fields = split(line, '\t')
            all_strains = String[String(fields[i]) for i in 10:length(fields)]
            break
        end
    end

    debug_log("VCF header: ", length(all_strains), " samples, ",
              length(chrom_rank), " contigs")
    (all_strains, chrom_rank, info_headers)
end

"""
    open_vcf_peeked(path) -> (PeekedFile, all_strains, chrom_rank, info_headers)

Opens a bgzip-compressed VCF via subprocess, parses its header, returns
a PeekedFile positioned at the first data line.
"""
function open_vcf_peeked(path::String)
    io = open(`bgzip -d -c $path`)
    (all_strains, chrom_rank, info_headers) = parse_vcf_header(io)
    pf = PeekedFile(io, "", false)
    advance!(pf)
    (pf, all_strains, chrom_rank, info_headers)
end
```

- [ ] **Step 4: Delete `peek_end_key` (lines 833-849)**

Delete the entire `peek_end_key` function and its docstring.

- [ ] **Step 5: Update section comment above `peek_sort_key`**

Change:
```julia
# Sorted-merge helpers: sort keys over VCF/GVCF lines
```
to:
```julia
# Sorted-merge helpers: sort keys over VCF lines
```

- [ ] **Step 6: Verify syntax**

```bash
docker run --rm -v $(pwd):/work -w /work veupathdb/dnaseqanalysis:1.0.0 \
  julia --compile=min -e 'include("bin/processSequenceVariations.jl")' 2>&1 | head -20
```

Expected: no parse errors.

- [ ] **Step 7: Commit**

```bash
git add bin/processSequenceVariations.jl
git commit -m "refactor: rename GVCFRecord->VCFRecord, remove is_ref_block/end_pos, remove peek_end_key"
```

---

### Task 9: Replace `prev_coverage_span` in `build_variations_from_record`

**Files:**
- Modify: `bin/processSequenceVariations.jl` — `build_variations_from_record` (lines ~739-805)

- [ ] **Step 1: Replace the full function**

```julia
"""
    build_variations_from_record(record, all_strains, undone_strains, chrom_coverage)
        -> Vector{Variation}

Builds per-strain Variation records from a VCF variant record.
For missing GTs, synthesizes a reference call if coverage.tsv shows the position covered.
"""
function build_variations_from_record(
    record::VCFRecord,
    all_strains::Vector{String},
    undone_strains::Set{String},
    chrom_coverage::Dict{String, Vector{Tuple{Int, Int, Float64}}}
)::Vector{Variation}
    variations = Variation[]

    for (i, strain) in enumerate(all_strains)
        strain in undone_strains && continue
        i > length(record.sample_data) && continue

        fmt = parse_format_field(record.format_keys, record.sample_data[i])

        gt = get(fmt, "GT", "")
        if isempty(gt) || gt == "." || gt == "./." || gt == ".|."
            # No call: synthesize a reference Variation if position is covered
            (covered, dp) = get_coverage(chrom_coverage, strain, record.pos - 1)
            if covered
                v = Variation()
                v.sequence_source_id = record.chrom
                v.location           = record.pos
                v.strain             = strain
                v.reference          = record.ref
                v.base               = record.ref
                v.coverage           = string(dp)
                v.percent            = "100"
                v.quality            = "."
                v.pvalue             = "."
                v.snp_source_id      = "NGS_SNP.$(record.chrom).$(record.pos)"
                v.matches_reference  = 1
                push!(variations, v)
            end
            continue
        end

        dp_str = get(fmt, "DP", "0")
        dp     = isempty(dp_str) || dp_str == "." ? 0 : parse(Int, dp_str)

        base = gt_to_base(gt, record.ref, record.alts)
        isempty(base) && continue

        aidx = gt_allele_idx(gt)
        pct  = compute_percent(fmt, aidx)
        gq   = get(fmt, "GQ", "0")

        v = Variation()
        v.sequence_source_id = record.chrom
        v.location           = record.pos
        v.strain             = strain
        v.reference          = record.ref
        v.base               = base
        v.coverage           = string(dp)
        v.percent            = pct
        v.quality            = gq
        v.pvalue             = "."
        v.snp_source_id      = "NGS_SNP.$(record.chrom).$(record.pos)"
        v.matches_reference  = (base == record.ref) ? 1 : 0

        push!(variations, v)
    end

    variations
end
```

Note: `record.pos - 1` converts 1-based VCF position to 0-based BED coordinate for coverage lookup.

- [ ] **Step 2: Verify syntax**

```bash
docker run --rm -v $(pwd):/work -w /work veupathdb/dnaseqanalysis:1.0.0 \
  julia --compile=min -e 'include("bin/processSequenceVariations.jl")' 2>&1 | head -20
```

Expected: no parse errors.

- [ ] **Step 3: Commit**

```bash
git add bin/processSequenceVariations.jl
git commit -m "refactor: replace prev_coverage_span with chrom_coverage in build_variations_from_record"
```

---

### Task 10: Replace `prev_coverage_span` in `handle_variant_record!` and update its signature

**Files:**
- Modify: `bin/processSequenceVariations.jl` — `handle_variant_record!` (lines ~1439-1585)

- [ ] **Step 1: Update the function signature (lines 1439-1447)**

Replace the docstring and signature:

```julia
"""
    handle_variant_record!(record, cache_entries, ctx, writers, transcript_cache, all_strains, chrom_coverage) -> Bool

Processes one variant VCF record end-to-end. Returns true if output was written.
"""
function handle_variant_record!(
    record::VCFRecord,
    cache_entries::Dict{Tuple{String,String},CacheEntry},
    ctx::ProcessingContext,
    writers::OutputWriters,
    transcript_cache::TranscriptSequenceCache,
    all_strains::Vector{String},
    chrom_coverage::Dict{String, Vector{Tuple{Int, Int, Float64}}}
)::Bool
```

- [ ] **Step 2: Update the `build_variations_from_record` call (line ~1452)**

Change:
```julia
    variations = build_variations_from_record(record, all_strains, ctx.undone_strains, prev_coverage_span)
```
to:
```julia
    variations = build_variations_from_record(record, all_strains, ctx.undone_strains, chrom_coverage)
```

- [ ] **Step 3: Replace the `modified_sample_data` span block (lines 1513-1531)**

Replace from `# Build modified sample data...` through the closing `end` of the per-strain loop with:

```julia
    # Build modified sample data: fill in GT=0 and DP for samples that are covered
    # at this position but were left as missing GT by bcftools merge.
    gt_idx = findfirst(==("GT"), record.format_keys)
    dp_idx = findfirst(==("DP"), record.format_keys)
    modified_sample_data = copy(record.sample_data)
    for (i, strain) in enumerate(all_strains)
        i > length(modified_sample_data) && continue
        fmt = parse_format_field(record.format_keys, modified_sample_data[i])
        gt = get(fmt, "GT", "")
        (isempty(gt) || gt == "." || gt == "./." || gt == ".|.") || continue
        (covered, dp) = get_coverage(chrom_coverage, strain, record.pos - 1)
        covered || continue
        fields = fill(".", length(record.format_keys))
        !isnothing(gt_idx) && (fields[gt_idx] = "0")
        !isnothing(dp_idx) && (fields[dp_idx] = string(round(Int, dp)))
        modified_sample_data[i] = join(fields, ":")
    end
```

- [ ] **Step 4: Verify syntax**

```bash
docker run --rm -v $(pwd):/work -w /work veupathdb/dnaseqanalysis:1.0.0 \
  julia --compile=min -e 'include("bin/processSequenceVariations.jl")' 2>&1 | head -20
```

Expected: no parse errors.

- [ ] **Step 5: Commit**

```bash
git add bin/processSequenceVariations.jl
git commit -m "refactor: replace prev_coverage_span with chrom_coverage in handle_variant_record!"
```

---

### Task 11: Refactor `main()` in `processSequenceVariations.jl`

**Files:**
- Modify: `bin/processSequenceVariations.jl` — `main()` function (lines 1591-1693)

- [ ] **Step 1: Update the file header comment (lines 1-8)**

```julia
#!/usr/bin/env julia

# processSequenceVariations.jl
# Reads a merged multi-sample FreeBayes VCF and a coordinate-sorted VCF cache file,
# streams them concurrently in a sorted merge, annotates coding variants via SQLite
# transcript/indel databases, and writes four output files:
#   cache.vcf (CANN-annotated VCF cache), snpFeature.dat, allele.dat, product.dat
# Coverage information is read from a coverage.tsv produced by mergeCoverageBeds.
```

- [ ] **Step 2: Replace the full `main()` function (lines 1591-1693)**

```julia
function main()
    args = parse_args(ARGS)
    global DEBUG = haskey(args, "debug")
    debug_log("Debug mode enabled")

    # Open VCF and parse header
    debug_log("Opening VCF: ", args["vcf_file"])
    (vcf_pf, all_strains, chrom_rank, info_headers) = open_vcf_peeked(args["vcf_file"])
    debug_log("VCF: ", length(all_strains), " strains")

    # Open VCF cache (may be absent/empty on first run)
    debug_log("Opening cache: ", args["cache_file"])
    cache_pf = open_cache_peeked(args["cache_file"])

    # Open coverage file
    debug_log("Opening coverage: ", args["coverage_file"])
    coverage_fh = open_coverage_file(args["coverage_file"])

    # Initialize processing context
    ctx = initialize_processing_context(args, all_strains)
    debug_log("Context: ", length(ctx.all_strains), " strains, ",
              length(ctx.cds_intervals), " CDS intervals")

    # Open output writers and write VCF cache header
    writers = open_output_writers(args["output_vcf"])
    write_vcf_cache_header(writers.vcf_cache_fh, ctx.all_strains, info_headers)

    transcript_cache = TranscriptSequenceCache(Dict{String, Dict{String,String}}())

    chrom_coverage = Dict{String, Vector{Tuple{Int, Int, Float64}}}()
    current_chrom  = ""

    n_processed = 0

    while !vcf_pf.exhausted
        vcf_start = peek_sort_key(vcf_pf.line, chrom_rank)
        cache_key = cache_pf.exhausted ? (typemax(Int), typemax(Int)) :
                                          peek_sort_key(cache_pf.line, chrom_rank)

        # Drain cache entries that precede the current VCF record start
        # (positions that were variant in a prior run but are now absent)
        if !cache_pf.exhausted && cache_key < vcf_start
            advance!(cache_pf)
            continue
        end

        # Parse and advance VCF
        record = parse_vcf_record(vcf_pf.line, length(all_strains))
        advance!(vcf_pf)

        # Load coverage intervals when the chromosome changes
        if record.chrom != current_chrom
            current_chrom = record.chrom
            load_chrom_coverage!(coverage_fh, current_chrom, chrom_rank, chrom_coverage)
        end

        # Collect all cache entries at this (chrom, pos)
        cache_entries = Dict{Tuple{String,String},CacheEntry}()
        while !cache_pf.exhausted
            ck = peek_sort_key(cache_pf.line, chrom_rank)
            ck != vcf_start && break
            parsed = parse_cache_vcf_record(cache_pf.line)
            if !isnothing(parsed)
                (_, _, ref, alt, cann_str) = parsed
                cache_entries[(ref, alt)] = CacheEntry(cann_str)
            end
            advance!(cache_pf)
        end

        if handle_variant_record!(record, cache_entries, ctx, writers, transcript_cache, all_strains, chrom_coverage)
            n_processed += 1
            if n_processed % 1000 == 0
                @info "Processed $n_processed variant positions"
            end
        end
    end

    debug_log("Processing complete. Total positions processed: ", n_processed)

    close_peeked(vcf_pf)
    close_peeked(cache_pf)
    close(coverage_fh.fh)
    close_output_writers(writers)
    close_processing_context(ctx)

    debug_log("Done!")
end
```

- [ ] **Step 3: Verify syntax and confirm no stale symbols remain**

```bash
docker run --rm -v $(pwd):/work -w /work veupathdb/dnaseqanalysis:1.0.0 \
  julia --compile=min -e 'include("bin/processSequenceVariations.jl")' 2>&1 | head -30
```

```bash
grep -n "prev_coverage_span\|is_ref_block\|peek_end_key\|GVCFRecord\|gvcf_pf\|open_gvcf\|parse_gvcf" \
  bin/processSequenceVariations.jl
```

Expected: Docker command shows no parse errors; grep returns no matches.

- [ ] **Step 4: Commit**

```bash
git add bin/processSequenceVariations.jl
git commit -m "feat: refactor main() to use coverage.tsv per-chrom streaming; remove REF block loop"
```

---

### Task 12: End-to-end stub verification

**Files:** none (verification only)

- [ ] **Step 1: Stub run `processSingleExperiment`**

```bash
nextflow run main.nf -profile processSingleExperiment -stub
```

Expected: all stages complete. `makeCoverageBed` present; `mergeCoverageBeds` absent.

- [ ] **Step 2: Stub run `mergeExperiments`**

```bash
nextflow run main.nf -entry mergeExperiments -profile mergeExperiments -stub
```

Expected: stages `mergeVcfs`, `mergeCoverageBeds`, `makeGenomicIndelDb`, `makeCodingData`, `processSeqVars`, `snpEff` all present. `mergeGvcfs` absent.

- [ ] **Step 3: Commit any remaining fixes**

```bash
git add -p
git commit -m "fix: resolve stub run issues"
```
