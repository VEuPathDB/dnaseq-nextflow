#!/usr/bin/env julia

# Characterization tests for the in-memory GenomicIndelIndex in makeCodingData.jl,
# which replaced two per-exon SQL queries (fasta_offset's SUM(shift) and
# project_indels_to_cds's range scan) — the makeCodingData bottleneck. The
# in-memory lookups MUST match the original SQL exactly; this builds a small
# genomic_indels table, runs the original queries directly as the oracle, and
# asserts agreement across many positions and exon ranges.
#
# Run with: julia testing/t/genomicIndelIndex.jl

using Test
using SQLite
using SQLite.DBInterface: execute

include(joinpath(@__DIR__, "../../bin/GtfUtils.jl"))
include(joinpath(@__DIR__, "../../bin/makeCodingData.jl"))

# --- oracles: the exact SQL the in-memory index replaced --------------------

function sql_fasta_offset(db::SQLite.DB, strain::String, seq_id::String, before_pos::Int)
    row = first(execute(db,
        "SELECT COALESCE(SUM(shift), 0) FROM genomic_indels
         WHERE strain = ? AND sequence_id = ? AND position < ?",
        [strain, seq_id, before_pos]))
    Int(row[1])
end

function sql_project_range(db::SQLite.DB, strain::String, seq_id::String, lo::Int, hi::Int)
    [(Int(r[1]), Int(r[2])) for r in execute(db,
        "SELECT position, shift FROM genomic_indels
         WHERE strain = ? AND sequence_id = ? AND position >= ? AND position <= ?
         ORDER BY position",
        [strain, seq_id, lo, hi])]
end

function make_genomic_indel_db(rows::Vector)
    db = SQLite.DB()
    execute(db, """CREATE TABLE genomic_indels (
        strain TEXT, sequence_id TEXT, position INTEGER, shift INTEGER)""")
    for (strain, seq_id, pos, shift) in rows
        execute(db, "INSERT INTO genomic_indels VALUES (?,?,?,?)", [strain, seq_id, pos, shift])
    end
    db
end

# Deliberately exercise: multiple sequences, out-of-order insertion, duplicate
# positions (two events at the same position), positive & negative shifts, and a
# second strain that must not leak into the first strain's index.
const ROWS = [
    ("strainA", "chr1", 100,  1),
    ("strainA", "chr1", 50,   3),
    ("strainA", "chr1", 100,  2),   # duplicate position
    ("strainA", "chr1", 200, -1),
    ("strainA", "chr2", 10,   6),
    ("strainA", "chr2", 40,  -4),
    ("strainB", "chr1", 75,  -3),   # different strain — must be isolated
]

@testset "fasta_offset matches SQL oracle at every position" begin
    db  = make_genomic_indel_db(ROWS)
    idx = load_strain_genomic_indels(db, "strainA")
    for seq_id in ("chr1", "chr2", "chrX")   # chrX absent → always 0
        for pos in 0:250
            @test fasta_offset(idx, seq_id, pos) == sql_fasta_offset(db, "strainA", seq_id, pos)
        end
    end
    close(db)
end

@testset "project range matches SQL oracle over many windows" begin
    db  = make_genomic_indel_db(ROWS)
    idx = load_strain_genomic_indels(db, "strainA")
    positions, shifts = idx.byseq["chr1"][1], idx.byseq["chr1"][2]
    for lo in 0:210, hi in lo:210
        # replicate the in-memory slice the way project_indels_to_cds does
        a = searchsortedfirst(positions, lo)
        b = searchsortedlast(positions, hi)
        mem = [(positions[i], shifts[i]) for i in a:b]
        @test mem == sql_project_range(db, "strainA", "chr1", lo, hi)
    end
    close(db)
end

@testset "load isolates the requested strain" begin
    db  = make_genomic_indel_db(ROWS)
    idxB = load_strain_genomic_indels(db, "strainB")
    # strainB has only chr1@75; strainA's chr1 rows must not appear
    @test idxB.byseq["chr1"][1] == [75]
    @test fasta_offset(idxB, "chr1", 100) == -3
    @test fasta_offset(idxB, "chr1", 75)  == 0     # 75 not < 75
    idxNone = load_strain_genomic_indels(db, "strainZ")
    @test isempty(idxNone.byseq)
    @test fasta_offset(idxNone, "chr1", 100) == 0
    close(db)
end
