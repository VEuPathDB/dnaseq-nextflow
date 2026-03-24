# Plan: Refactor processSequenceVariations.jl — GVCF + VCF cache

## Context

The `processSeqVars` Nextflow process now passes a merged multi-sample GVCF (per-sample FreeBayes `--gvcf` output, merged by `bcftools merge`) instead of a custom 10-col SNP TSV. The per-strain coverage directory is eliminated — coverage (DP) is available directly from GVCF FORMAT fields. The old 20-col TSV annotation cache is replaced by a coordinate-sorted VCF file with a custom `CANN` INFO field (see CANN spec in memory). Both the GVCF and the VCF cache are coordinate-sorted, so they are streamed concurrently in a sorted-merge, analogous to the old SNP TSV + cache TSV merge. The VCF cache is written incrementally during the merge (not accumulated in memory). First run uses an empty cache file.

## Files to Modify

- `bin/processSequenceVariations.jl` — all changes here
- `nextflow.config` — add `minCoverage` to `mergeExperiments` profile; rename `vcfCacheFile = "cache.txt"` → `"cache.vcf"`
- `modules/mergeExperiments.nf` — add `--min_coverage $params.minCoverage` to `processSeqVars` script invocation

## Dead Code to Remove

| Symbol | Why removed |
|---|---|
| `CoverageReader` struct + all methods | Replaced by DP from GVCF FORMAT |
| `open_coverage_files` | Directory scan gone |
| `parse_cache_line` | 20-col TSV cache gone |
| `parse_snp_line` | 10-col SNP TSV gone |
| `write_cache_entries` | TSV cache writer gone; VCF cache written via new functions |
| `finalize_output_files` | SNP/undone_strains file truncation gone |
| `collect_variations_at_position` | Replaced by GVCF-aware version |
| `determine_next_position` | Replaced by GVCF/cache-aware version |
| `process_all_positions!` | Replaced by new main loop |
| `fill_coverage_gaps!` | Implicit in GVCF multi-sample column iteration |

## Kept / Adapted

- **`PeekedFile` struct** — kept as-is. Repurposed: one `PeekedFile` for the GVCF data stream, one for the VCF cache stream (instead of SNP TSV + cache TSV).
- **`open_peeked`, `advance!`, `close_peeked`** — kept unchanged.
- **`peek_sort_key`** — adapted for VCF lines: returns `(chrom_rank, pos)` from CHROM and POS columns.
- **`peek_end_key`** — **new helper** needed for REF block span awareness (see below).

## New/Updated Structs

### `GVCFRecord` (new)
Parsed representation of one GVCF data line, used in `handle_variant_record!`:
```julia
struct GVCFRecord
    chrom::String
    pos::Int
    ref::String
    alts::Vector{String}        # ["<*>"] for REF blocks
    is_ref_block::Bool
    format_keys::Vector{String}
    sample_data::Vector{String} # raw per-sample strings, order = all_strains
end
```

### `ProcessingContext` (updated)
- Remove `coverage_readers`
- Add `min_coverage::Int`
- `all_strains` passed in from GVCF header (not derived from coverage dir)

### `OutputWriters` (updated)
- Replace `cache_fh::IOStream` with `vcf_cache_fh::IOStream` — writes new VCF cache lines incrementally

## New Functions

### GVCF I/O
- **`open_gvcf(path) -> IO`** — `open(\`bgzip -d -c $path\`)` subprocess pipe (bgzip available in container)
- **`parse_gvcf_header(io) -> (Vector{String}, Dict{String,Int})`** — read `##` meta lines; build `chrom_rank` dict from `##contig` lines in order; on `#CHROM` line extract sample names (cols 10+); leaves io at first data line; returns `(all_strains, chrom_rank)`
- **`parse_gvcf_record(line) -> GVCFRecord`** — split on `\t`; detect REF block (`ALT == "<*>"` or symbolic); parse INFO/END for `end_pos` on REF blocks (else `end_pos = pos`); parse FORMAT column; store raw sample strings
- **`parse_format_field(format_keys, sample_str) -> Dict{String,String}`** — zip FORMAT keys with sample colon-fields; graceful if sample has fewer values than keys
- **`peek_end_key(pf::PeekedFile, chrom_rank) -> (Int, Int)`** — **new**: for REF blocks returns `(rank, END)` from INFO field; for variant records returns `(rank, POS)`. Used in sorted-merge to know how far a REF block spans without fully parsing it. Implementation: split peeked line on `\t`, check if ALT is symbolic, if so regex `END=(\d+)` from INFO field.

### VCF Cache I/O
- **`open_cache(path) -> IO`** — plain `open(path, "r")`; if file absent/empty return a pre-exhausted `PeekedFile`
- **`parse_cache_vcf_record(line) -> (chrom, pos, ref, alt, cann_str)`** — parse one VCF cache line; extract CANN value from INFO field (`strip("CANN=", fields[8])`)
- **`write_vcf_cache_header(fh)`** — writes `##fileformat=VCFv4.2`, `##INFO=<ID=CANN,...>`, `#CHROM...`
- **`write_vcf_cache_entry(fh, chrom, pos, ref, alt, cann_str)`** — one tab-sep line

### CANN Construction
- **`build_cann_string(outcomes::Vector) -> String`** — dedup `(codon, alt_aa, effect, transcript_id, pos_in_cds, pos_in_codon)` tuples → assign k0/k1/..., comma-join. `.` for absent fields. Non-coding → `"."`.
- **`classify_effect(...) -> String`** — per CANN spec table: `missense`, `synonymous`, `nonsense`, `frameshift`, `upstream_frameshift`, `inframe_insertion`, `inframe_deletion`, compound with `;`
- **`decode_cann_to_annotation(cann_str) -> PositionAnnotation or nothing`** — parse k0 entry to reconstruct `is_coding`, `transcript_id`, `pos_in_cds`, `pos_in_codon`, `ref_codon`; `nothing` for `"."`

### Variation Construction from GVCF
- **`gt_to_base(gt, ref, alts) -> String`** — haploid `"0"/"1"`, diploid `"0/0"/"0/1"/"1/1"`, phased `|`. Het SNP → IUPAC code via existing table.
- **`compute_percent(fmt_dict, allele_idx) -> String`** — `AO[allele_idx] / (RO + AO[allele_idx]) * 100`; `"0.0"` for ref or missing
- **`build_variations_from_record(record, all_strains, undone_strains, min_coverage) -> Vector{Variation}`** — iterate all sample columns; skip `"./."` and DP < min_coverage; GT=`0/0` → ref-call Variation; GT with alt → alt Variation. Synthesizes `snp_source_id = "NGS_SNP.$(chrom).$(pos)"`, `pvalue = "."`, `quality = get(fmt, "GQ", "0")`

### Core Record Handler
- **`handle_variant_record!(gvcf_record, cann_str_or_nothing, ctx, writers, transcript_cache) -> Bool`** — replaces `process_single_position!`:
  1. `build_variations_from_record` → variations
  2. If `cann_str_or_nothing != nothing` → `decode_cann_to_annotation` (skip `annotate_position`)
     else → `annotate_position` (full CDS lookup)
  3. `annotate_variations!` — always runs (per-strain indel DB queries not cached)
  4. Set `downstream_of_frameshift` on each Variation (moved from old `write_cache_entries`)
  5. `build_reference_variation`
  6. `has_variation` check → return false if all ref
  7. `build_cann_string` per ALT allele
  8. `write_vcf_cache_entry` for each (chrom, pos, ref, alt) → incremental cache write
  9. `write_snp_feature`, `write_allele_and_product_files`
  10. return true

## `annotate_variations!` Change

Remove the guard `v.position_in_codon != 0` (~line 1011) that skipped already-annotated Variations. Cache hits are now handled at the position level (skip `annotate_position`), not per-Variation. All Variations arrive with `position_in_codon = 0`.

## Updated `main()` Flow — Sorted Concurrent Stream

```
1.  parse_args(ARGS)
2.  gvcf_io = open_gvcf(args["vcf_file"])
3.  (all_strains, chrom_rank) = parse_gvcf_header(gvcf_io)
    # io now positioned at first data line
4.  gvcf_pf  = PeekedFile(gvcf_io, ...)   # peek first data line
5.  cache_pf = open_peeked(args["cache_file"])   # empty PeekedFile if file absent

6.  ctx = initialize_processing_context(args, all_strains)

7.  (writers, temp_cache_file) = open_output_writers(args["cache_file"])
    write_vcf_cache_header(writers.vcf_cache_fh)

8.  transcript_cache = TranscriptSequenceCache("", Dict())

9.  Main loop — span-aware sorted merge of gvcf_pf and cache_pf:
    while !gvcf_pf.exhausted:
        gvcf_start = peek_sort_key(gvcf_pf, chrom_rank)    # (rank, POS)
        gvcf_end   = peek_end_key(gvcf_pf, chrom_rank)     # (rank, END) for REF blocks; (rank, POS) for variants
        cache_key  = cache_pf.exhausted ? (∞,∞) : peek_sort_key(cache_pf, chrom_rank)

        # Drain cache entries that precede the current GVCF record's START
        # (positions no longer in GVCF from a prior run)
        if cache_key < gvcf_start:
            advance!(cache_pf)
            continue

        record = parse_gvcf_record(gvcf_pf.line)
        advance!(gvcf_pf)

        if record.is_ref_block:
            # Explicitly drain ALL cache entries within the REF block span [POS..END]
            # These positions are now reference in this run — drop from new cache
            while !cache_pf.exhausted && peek_sort_key(cache_pf, chrom_rank) <= gvcf_end:
                advance!(cache_pf)
            continue

        # Variant record — check for cache hit
        if cache_key == gvcf_start:
            (_, _, _, _, cann_str) = parse_cache_vcf_record(cache_pf.line)
            advance!(cache_pf)
            handle_variant_record!(record, cann_str, ctx, writers, transcript_cache)
        else:
            handle_variant_record!(record, nothing, ctx, writers, transcript_cache)

10. close(gvcf_io); close_peeked(cache_pf)
11. close_output_writers(writers)
12. close_processing_context(ctx)
13. mv(temp_cache_file, args["cache_file"])
```

**REF block span handling**: `peek_end_key` exposes the INFO/END coordinate of REF blocks without a full parse. When a REF block is encountered, the inner `while` loop explicitly drains all cache entries within `[POS, END]`. This is correct and explicit — those positions were variant in a prior run but are now reference-covered; they must not appear in the new cache.

**Sort key comparison**: `peek_sort_key` returns `(chrom_rank, pos)`. `chrom_rank` is built from `##contig` lines in the GVCF header, ensuring GVCF and VCF cache use identical chromosome ordering (guaranteed since the VCF cache was produced by a prior run of this script on the same genome).

**Cache write is incremental**: `write_vcf_cache_entry` is called inside `handle_variant_record!` after CANN is computed. No in-memory accumulation.

**Multi-allelic same-position cache entries**: If the VCF cache has two entries at the same (chrom, pos) — one for each ALT of a prior multi-allelic call — the inner `cache_key == gvcf_start` loop must drain ALL cache entries at that position (not just the first). Collect them into a `Dict{(ref,alt) => cann_str}` and pass the whole dict into `handle_variant_record!`.

## Chromosome Sort Key

`parse_gvcf_header` reads `##contig=<ID=chr1,...>` lines in order to build:
```julia
chrom_rank = Dict{String,Int}("chr1"=>1, "chr2"=>2, ...)
```
`peek_sort_key(line, chrom_rank) -> (Int, Int)`: split line on `\t`, return `(chrom_rank[fields[1]], parse(Int, fields[2]))`. If chrom not in rank dict (shouldn't happen), fall back to hash. Both GVCF and VCF cache must use the same ordering — guaranteed since the cache was written by a prior run of this same script from the same GVCF chromosome order.

## Arg Changes

| Arg | Change |
|---|---|
| `--snp_file` | **Removed** |
| `--coverage_directory` | **Removed** |
| `--vcf_file` | Kept (already in Nextflow process) |
| `--min_coverage` | **Added** (default 1) |
| All others | Unchanged |

## Edge Cases

- **Haploid GT**: FreeBayes `--ploidy 1` produces `"0"/"1"`. `gt_to_base` handles single-value GT (no `/`).
- **Multi-allelic ALT**: FreeBayes can produce multi-allelic records; AO is comma-sep. GT allele index maps to `alts[idx]`. One cache entry per (chrom, pos, ref, alt) pair.
- **Missing FORMAT fields**: GQ absent in some FreeBayes records. `get(fmt, "GQ", "0")`.
- **Non-coding in cache**: `CANN=.` stored so re-runs skip CDS lookup for known non-coding positions.
- **Empty cache first run**: `open_peeked` on absent/empty file returns a pre-exhausted `PeekedFile`; loop falls entirely into the `gvcf_key < cache_key` branch.
- **`downstream_of_frameshift`**: Moved from old `write_cache_entries` to `handle_variant_record!`; must be set on each Variation before calling `write_allele_and_product_files`.
- **`vcfCacheFile` rename**: `"cache.txt"` → `"cache.vcf"` in config.
- **Multi-allelic cache key collision**: If a GVCF record has REF=A, ALT=T,G and the cache has entries for both (A,T) and (A,G), the sort key `(chrom, pos)` matches both. Handle by draining all cache entries at the same (chrom, pos) into a Dict keyed by (ref, alt) before calling `handle_variant_record!`.

## Unchanged (Do Not Touch)

`annotate_position`, `annotate_variations!` (minus cache guard), codon translation, IUPAC expansion, GTF parsing, `precompute_frameshifts`, `TranscriptSequenceCache`, `build_reference_variation`, `has_variation`, `write_snp_feature`, `write_allele_and_product_files`, `CDSInterval`, `TranscriptInfo`, `PositionAnnotation`, `Variation` structs.

## Implementation Order

1. GVCF parsing (`open_gvcf`, `parse_gvcf_header`, `parse_gvcf_record`, `parse_format_field`)
2. VCF cache I/O (`open_cache`, `parse_cache_vcf_record`, `write_vcf_cache_header`, `write_vcf_cache_entry`)
3. Variation construction (`gt_to_base`, `compute_percent`, `build_variations_from_record`)
4. CANN layer (`classify_effect`, `build_cann_string`, `decode_cann_to_annotation`)
5. Remove dead code; update `ProcessingContext`, `OutputWriters`, `initialize_processing_context`
6. Update `annotate_variations!` (remove cache guard)
7. Update `peek_sort_key` for VCF lines + chrom_rank; add `peek_end_key` for REF block span
8. `handle_variant_record!`
9. Update `main()` (span-aware sorted concurrent stream loop)

## Verification

```bash
# Stub run (confirms Nextflow wiring with no real data):
nextflow run main.nf -entry mergeExperiments -profile mergeExperiments -stub

# Full run with test data:
nextflow run main.nf -entry mergeExperiments -profile mergeExperiments

# Manual Julia smoke test:
# julia bin/processSequenceVariations.jl \
#   --vcf_file data/merge_setup/merged.vcf.gz \
#   --cache_file data/merge_setup/cache.vcf \
#   --undone_strains_file data/merge_setup/undoneStrains.txt \
#   --reference_strain pfal3D7 \
#   --transcript_db data/codingSequences.db \
#   --indel_db data/codingIndels.db \
#   --gtf_file data/ref.gtf \
#   --min_coverage 5
```
