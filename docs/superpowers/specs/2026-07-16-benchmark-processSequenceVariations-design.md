# Benchmark instrumentation for `processSequenceVariations.jl`

**Date:** 2026-07-16
**Status:** Approved design, pending implementation
**Scope:** Measure only. No behavior change; no hotspot fix in this change.

## Problem

`mergeExperiments` takes ~14 hours for 200 DNA-seq samples, and most of that
time is spent inside `bin/processSequenceVariations.jl`. Production runs for some
organisms process 1000–2000 samples, so whatever scales poorly with sample count
dominates. We need to know *where* the time goes before optimizing.

A complication: the only convenient local dataset has 4 samples. Wall-clock on 4
samples is small, so raw timing alone will not surface an O(samples) hotspot. The
instrumentation must therefore also expose **per-position call counts**, which
reveal per-sample scaling regardless of absolute time.

## Prime suspects (from code reading)

- `get_indel_shift` (`bin/processSequenceVariations.jl`) issues one SQLite query
  **per variation per coding transcript**, inside `annotate_variations!`.
- Reference synthesis in `build_variations_from_record` creates a `Variation` for
  **every covered sample** at **every no-call position**.

Together these make coding loci roughly O(samples) or worse in SQL round-trips.
This design measures them precisely rather than assuming.

## Mechanism

A global `BENCHMARK` flag parallel to the existing `DEBUG` flag, set from a new
`--benchmark` CLI argument. A `@bench` macro wraps an expression:

```julia
macro bench(name, expr)
    quote
        if BENCHMARK
            local t0 = time_ns()
            local v  = $(esc(expr))
            _bench_record($(esc(name)), time_ns() - t0)
            v
        else
            $(esc(expr))
        end
    end
end
```

When `BENCHMARK` is false the macro expands to just `expr` — zero overhead on
normal runs. `_bench_record(name, ns)` accumulates `(total_ns, count)` into a
global `Dict{String,Tuple{Int,Int}}` (ns as Int, count as Int).
`bench_count!(name, n=1)` bumps a pure counter (no timing) in the same store.

## Instrumentation points

Setup (once):
- `parse_gtf`
- `precompute_frameshifts`
- VCF open + header parse

Per-position phases (inside `handle_variant_record!` and the main loop):
- `build_variations`
- `annotate` (covers both the cache path and the fresh GTF lookup)
- `annotate_variations`
- `write_allele`
- `write_transcript_product`
- `write_snp_feature`
- `cann_collect` (CANN entry collection + `assign_cann_keys`)
- `write_hsss`
- `vcf_output_write`

Hot SQL calls, wrapped individually for time **and** count:
- `sql_get_indel_shift`
- `sql_load_transcript`

Scale counters (count only):
- `positions_processed`
- `records_parsed`
- `variations_built` (total `Variation` objects created — reveals synthesis blowup)

## Output

Enabled by `--benchmark`. At the end of `main()`, print a formatted table to
**stderr** (Nextflow captures stderr in `.command.err`, so it is visible in
production runs too). Columns:

| bucket | total time | % of run | calls | avg µs/call | calls/position |

Rows sorted by total time descending. A header line reports total wall-clock and
`positions_processed` so the per-position column is interpretable.

The `calls/position` column is the key cross-scale signal: e.g. if
`sql_get_indel_shift` shows 8 calls/position on 4 samples (~2 per sample per
position), a 2000-sample run implies ~4000 SQL round-trips per coding locus —
legible even though the local wall-clock is small.

## Non-goals

- No optimization / hotspot fix in this change (separate, evidence-backed change).
- No new output artifacts, so the Nextflow process definition is untouched.
- No new package dependencies (`time_ns` is in `Base`).

## Testing

- Zero-overhead-when-off is structural (macro expands to the bare expression);
  confirmed by inspection.
- Run `--benchmark` against the 4-sample local dataset inside the
  `dnaseqanalysis` container and confirm the report prints with sane numbers and
  the phase percentages sum to ~100%.
- Existing Julia/Python test suites must still pass unchanged.

## Findings (20-sample run, 122,530 positions)

Benchmark on a real 20-sample Plasmodium dataset took **1661.84 s (~28 min)** and
was ~95% SQLite:

| bucket | time | note |
|---|---|---|
| `sql_load_transcript` | 931 s (56%) | only 4034 calls but **231 ms each** |
| `sql_get_indel_shift` | 627 s (38%) | **778,513 calls** (6.35/pos, linear in samples) at ~800 µs each |

Two distinct pathologies, confirmed against the actual DBs:

1. `coding_sequences` has `PRIMARY KEY (strain, transcript_id)`. The query filters
   `WHERE transcript_id = ?` only, so the autoindex (leading with `strain`) cannot
   serve it — every call full-scans the 364 MB table.
2. The `indels` table is tiny (5341 rows) and already indexed; the cost was purely
   778k statement round-trips, and the call count grows with sample count.

The output writers originally suspected were <2% combined.

## Fixes implemented (follow-up to the measure-only scope)

1. **`bin/makeCodingData.jl`** — `finalize_cds_db` adds
   `CREATE INDEX idx_coding_sequences_transcript ON coding_sequences(transcript_id)`
   at DB-build time. The read-only consumer is not mutated.
2. **`bin/processSequenceVariations.jl`** — `precompute_indel_shifts` scans the
   `indels` table once into a per-`(strain, transcript_id)` prefix-sum table;
   `lookup_indel_shift` answers `position < p` by binary search. `get_indel_shift`
   (the per-variation SQL query) is removed.
3. **`testing/t/indelShiftLookup.jl`** — characterization test asserting
   `lookup_indel_shift` equals the original SQL at every position (1755 checks + edge cases).

## Verification

Re-ran the fixed code against the same 20-sample inputs (with the index added to a
copy of the DB):

- **1661.84 s → 84.07 s (~20× faster).** `sql_load_transcript` 931 s → 1.9 s;
  `sql_get_indel_shift` 627 s → eliminated (0.27 s one-time precompute).
- **All outputs byte-identical** (variationFeature/snpFeature.dat,
  transcript_product.dat, allele.dat, sample.dat, output.vcf, all four HSSS
  binary dirs). Full Julia test suite passes.

Both costs were O(samples), so the projected effect at 1000–2000 samples is far
larger than 20× (the index fix removes a per-call cost that itself grew with
table size).
