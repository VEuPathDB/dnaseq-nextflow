#!/usr/bin/env julia

# Characterization tests for precompute_indel_shifts / lookup_indel_shift in
# processSequenceVariations.jl. These replaced a per-variation SQL aggregate
# (SUM(shift_amount) WHERE transcript_id=? AND strain=? AND position < ?) with a
# precomputed prefix-sum + binary search. The lookup MUST match the old SQL
# semantics exactly — this test builds a small indels table, runs the original
# query directly, and asserts the in-memory lookup agrees at every position.
#
# Run with: julia testing/t/indelShiftLookup.jl

using Test
using SQLite
using SQLite.DBInterface: execute

include(joinpath(@__DIR__, "../../bin/processSequenceVariations.jl"))

# The exact SQL the precompute replaced — kept here as the oracle.
function sql_indel_shift(db::SQLite.DB, tid::String, strain::String, position::Int)
    row = first(execute(db,
        "SELECT COALESCE(SUM(shift_amount), 0) FROM indels WHERE transcript_id = ? AND strain = ? AND position < ?",
        [tid, strain, position]))
    row[1]::Int
end

function make_indels_db(rows::Vector)
    db = SQLite.DB()
    execute(db, """
        CREATE TABLE indels (
          strain        TEXT    NOT NULL,
          transcript_id TEXT    NOT NULL,
          position      INTEGER NOT NULL,
          shift_amount  INTEGER NOT NULL
        )
    """)
    for (strain, tid, pos, shift) in rows
        execute(db, "INSERT INTO indels VALUES (?,?,?,?)", [strain, tid, pos, shift])
    end
    db
end

# Rows deliberately exercise: multiple strains & transcripts, out-of-order
# insertion, duplicate positions (two events at the same position), and both
# positive (insertion) and negative (deletion) shifts.
const ROWS = [
    ("strainA", "T1", 100,  1),
    ("strainA", "T1", 50,   3),
    ("strainA", "T1", 100,  2),   # duplicate position with T1@100 above
    ("strainA", "T1", 200, -1),
    ("strainA", "T2", 10,   6),
    ("strainB", "T1", 75,  -3),
    ("strainB", "T1", 300,  9),
]

@testset "lookup_indel_shift matches SQL oracle at every position" begin
    db     = make_indels_db(ROWS)
    shifts = precompute_indel_shifts(db)

    # Every (strain, transcript) pair present, plus a pair with no indels.
    pairs = [("strainA", "T1"), ("strainA", "T2"), ("strainB", "T1"),
             ("strainA", "T3"), ("strainC", "T1")]

    for (strain, tid) in pairs
        for position in 0:350
            expected = sql_indel_shift(db, tid, strain, position)
            actual   = lookup_indel_shift(shifts, tid, strain, position)
            @test actual == expected
        end
    end
    close(db)
end

@testset "lookup_indel_shift specific values" begin
    db     = make_indels_db(ROWS)
    shifts = precompute_indel_shifts(db)

    # strainA/T1 events: pos 50(+3), 100(+1), 100(+2), 200(-1)
    @test lookup_indel_shift(shifts, "T1", "strainA", 50)  == 0    # nothing strictly before 50
    @test lookup_indel_shift(shifts, "T1", "strainA", 51)  == 3    # only pos 50
    @test lookup_indel_shift(shifts, "T1", "strainA", 100) == 3    # 100 not < 100
    @test lookup_indel_shift(shifts, "T1", "strainA", 101) == 6    # 50 + both @100 (3+1+2)
    @test lookup_indel_shift(shifts, "T1", "strainA", 999) == 5    # all: 3+1+2-1
    @test lookup_indel_shift(shifts, "T1", "strainA", 0)   == 0

    # Absent pair → 0
    @test lookup_indel_shift(shifts, "T9", "strainZ", 100) == 0
    close(db)
end
