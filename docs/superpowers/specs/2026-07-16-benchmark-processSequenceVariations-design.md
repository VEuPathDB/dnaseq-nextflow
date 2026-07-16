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
