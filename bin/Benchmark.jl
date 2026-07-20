# Benchmark.jl — shared lightweight instrumentation for the pipeline's Julia
# scripts (included via `include("$(@__DIR__)/Benchmark.jl")`).
#
# Enabled by a `--benchmark` CLI flag; the including script sets the global
# BENCHMARK accordingly. Zero overhead when BENCHMARK is false: `@bench name expr`
# expands to the bare `expr`. When enabled, it accumulates (total_ns, count) per
# named bucket so an end-of-run report can attribute both wall-clock time AND
# call counts — the latter exposes work that scales with sample/strain count even
# on a tiny local dataset where absolute time is small.

using Printf

global BENCHMARK = false

# name -> [total_ns, count]
const BENCH_DATA = Dict{String, Vector{Int}}()

function _bench_record(name::AbstractString, ns::Integer)
    d = get!(BENCH_DATA, name, Int[0, 0])
    d[1] += ns
    d[2] += 1
    nothing
end

# Bump a pure counter (no timing) — for scale metrics like variations_built.
function bench_count!(name::AbstractString, n::Integer=1)
    BENCHMARK || return nothing
    d = get!(BENCH_DATA, name, Int[0, 0])
    d[2] += n
    nothing
end

macro bench(name, expr)
    quote
        if BENCHMARK
            local t0 = time_ns()
            local v  = $(esc(expr))
            _bench_record($(esc(name)), Int(time_ns() - t0))
            v
        else
            $(esc(expr))
        end
    end
end

"""
    print_benchmark_report(io, total_ns; summary, per_label, per_denom)

Print the accumulated buckets sorted by total time. `summary` is a free-text
header line (e.g. counts of the units processed). The last column normalizes each
bucket's call count by `per_denom` and is titled `per_label` (e.g. "calls/pos",
"calls/strain") — the signal that reveals per-unit scaling.
"""
function print_benchmark_report(io::IO, total_ns::Integer;
                                summary::AbstractString="",
                                per_label::AbstractString="calls/item",
                                per_denom::Integer=0)
    println(io, "")
    println(io, "==================== BENCHMARK REPORT ====================")
    @printf(io, "total wall-clock: %.2f s", total_ns / 1e9)
    isempty(summary) || print(io, "   ", summary)
    println(io)
    @printf(io, "%-26s %10s %6s %12s %10s %11s\n",
            "bucket", "total(s)", "%", "calls", "us/call", per_label)
    println(io, "-"^80)
    rows = sort(collect(BENCH_DATA); by = r -> -r[2][1])  # by total_ns desc
    for (name, d) in rows
        tns, cnt = d[1], d[2]
        pct     = total_ns > 0 ? tns / total_ns * 100 : 0.0
        us_call = cnt > 0 ? tns / cnt / 1e3 : 0.0
        per     = per_denom > 0 ? cnt / per_denom : 0.0
        @printf(io, "%-26s %10.3f %6.1f %12d %10.2f %11.3f\n",
                name, tns / 1e9, pct, cnt, us_call, per)
    end
    println(io, "==========================================================")
end
