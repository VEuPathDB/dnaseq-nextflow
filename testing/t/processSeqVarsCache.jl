#!/usr/bin/env julia

# Tests for TSV cache I/O in processSequenceVariations.jl
# Run with: julia testing/t/processSeqVarsCache.jl

using Test
using SQLite
using SQLite.DBInterface: execute

include(joinpath(@__DIR__, "../../bin/processSequenceVariations.jl"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function make_transcript_db(rows::Vector)
    db = SQLite.DB()
    execute(db, "CREATE TABLE coding_sequences (strain TEXT, transcript_id TEXT, sequence TEXT)")
    for (strain, tid, seq) in rows
        execute(db, "INSERT INTO coding_sequences VALUES (?,?,?)", [strain, tid, seq])
    end
    db
end

# ---------------------------------------------------------------------------
# parse_cache_tsv_record
# ---------------------------------------------------------------------------

@testset "parse_cache_tsv_record valid line" begin
    result = parse_cache_tsv_record("chr1\t100\tT1\t42")
    @test !isnothing(result)
    (chrom, pos, tid, pic) = result
    @test chrom == "chr1"
    @test pos   == 100
    @test tid   == "T1"
    @test pic   == 42
end

@testset "parse_cache_tsv_record comment line" begin
    @test isnothing(parse_cache_tsv_record("# header"))
end

@testset "parse_cache_tsv_record too few fields" begin
    @test isnothing(parse_cache_tsv_record("chr1\t100\tT1"))
end

@testset "parse_cache_tsv_record empty line" begin
    @test isnothing(parse_cache_tsv_record(""))
end

# ---------------------------------------------------------------------------
# open_cache_peeked (TSV)
# ---------------------------------------------------------------------------

@testset "open_cache_peeked absent file returns exhausted" begin
    pf = open_cache_peeked("/nonexistent/path/cache.tsv")
    @test pf.exhausted
end

@testset "open_cache_peeked empty file returns exhausted" begin
    tmp = tempname()
    touch(tmp)
    pf = open_cache_peeked(tmp)
    @test pf.exhausted
    rm(tmp)
end

@testset "open_cache_peeked skips comment header" begin
    tmp = tempname()
    open(tmp, "w") do f
        println(f, "# chrom\tpos\ttranscript_id\tpos_in_cds")
        println(f, "chr1\t50\tTX1\t10")
    end
    pf = open_cache_peeked(tmp)
    @test !pf.exhausted
    @test startswith(pf.line, "chr1")
    close_peeked(pf)
    rm(tmp)
end

# ---------------------------------------------------------------------------
# write_cache_entries
# ---------------------------------------------------------------------------

@testset "write_cache_entries coding annotations" begin
    buf = IOBuffer()
    annotations = [
        PositionAnnotation(1, "TX1", 0, 42, 0, "ATG", "M"),
        PositionAnnotation(1, "TX2", 0, 90, 0, "TTT", "F"),
    ]
    write_cache_entries(buf, "chr1", 100, annotations)
    lines = split(String(take!(buf)), '\n'; keepempty=false)
    @test length(lines) == 2
    @test lines[1] == "chr1\t100\tTX1\t42"
    @test lines[2] == "chr1\t100\tTX2\t90"
end

@testset "write_cache_entries skips non-coding" begin
    buf = IOBuffer()
    annotations = [PositionAnnotation(0, "", 0, 0, 0, "", "")]
    write_cache_entries(buf, "chr1", 100, annotations)
    @test isempty(String(take!(buf)))
end

# ---------------------------------------------------------------------------
# build_annotations_from_cache
# ---------------------------------------------------------------------------

@testset "build_annotations_from_cache single transcript" begin
    db = make_transcript_db([("reference", "TX1", "ATGTTTGGG")])
    tcache = TranscriptSequenceCache(Dict())
    ctx_ref = "reference"

    entries = [("TX1", 1)]
    anns = build_annotations_from_cache(entries, ctx_ref, db, tcache)
    @test length(anns) == 1
    @test anns[1].transcript_id == "TX1"
    @test anns[1].pos_in_cds    == 1
    @test anns[1].is_coding     == 1
    @test anns[1].ref_codon     == "ATG"
    @test anns[1].ref_product   == "M"
    @test anns[1].pos_in_codon_val == position_in_codon(1)
end

@testset "build_annotations_from_cache empty entries returns non-coding" begin
    db = make_transcript_db([])
    tcache = TranscriptSequenceCache(Dict())
    anns = build_annotations_from_cache(Tuple{String,Int}[], "reference", db, tcache)
    @test length(anns) == 1
    @test anns[1].is_coding == 0
end

@testset "build_annotations_from_cache multiple transcripts" begin
    db = make_transcript_db([
        ("reference", "TX1", "ATGTTTGGG"),
        ("reference", "TX2", "GGGCCCAAA"),
    ])
    tcache = TranscriptSequenceCache(Dict())
    entries = [("TX1", 1), ("TX2", 4)]
    anns = build_annotations_from_cache(entries, "reference", db, tcache)
    @test length(anns) == 2
    @test anns[1].transcript_id == "TX1"
    @test anns[2].transcript_id == "TX2"
    @test anns[2].pos_in_cds    == 4
    @test anns[2].ref_codon     == "CCC"
end

# ---------------------------------------------------------------------------
# Round-trip: write then read back
# ---------------------------------------------------------------------------

@testset "cache round-trip: write entries then parse them back" begin
    buf = IOBuffer()
    annotations = [
        PositionAnnotation(1, "TX1", 0, 1,  0, "ATG", "M"),
        PositionAnnotation(1, "TX1", 0, 10, 0, "TTT", "F"),
    ]
    write_cache_entries(buf, "chr1", 5,  [annotations[1]])
    write_cache_entries(buf, "chr1", 20, [annotations[2]])

    content = String(take!(buf))
    lines = split(content, '\n'; keepempty=false)
    @test length(lines) == 2

    r1 = parse_cache_tsv_record(lines[1])
    @test !isnothing(r1)
    @test r1 == ("chr1", 5, "TX1", 1)

    r2 = parse_cache_tsv_record(lines[2])
    @test !isnothing(r2)
    @test r2 == ("chr1", 20, "TX1", 10)
end
