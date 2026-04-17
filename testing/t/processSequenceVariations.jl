#!/usr/bin/env julia

using Test
using SQLite
using SQLite.DBInterface: execute

include(joinpath(@__DIR__, "../../bin/processSequenceVariations.jl"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function make_coding_indels_db(rows::Vector)
    db = SQLite.DB()
    execute(db, """CREATE TABLE indels (
        strain TEXT, transcript_id TEXT, position INTEGER,
        shift_amount INTEGER, zygosity TEXT, ref_allele TEXT, alt_allele TEXT)""")
    for (strain, tid, pos, shift, zyg, ref_a, alt_a) in rows
        execute(db, "INSERT INTO indels VALUES (?,?,?,?,?,?,?)",
                [strain, tid, pos, shift, zyg, ref_a, alt_a])
    end
    db
end

# ---------------------------------------------------------------------------
# lookup_cds_indel
# ---------------------------------------------------------------------------

@testset "lookup_cds_indel returns indel record" begin
    db = make_coding_indels_db([
        ("strainA", "T1", 10, -3, "hom", "ATGC", "A"),
    ])
    result = lookup_cds_indel(db, "T1", "strainA", 10)
    @test result !== nothing
    @test result == ("hom", "ATGC", "A")
end

@testset "lookup_cds_indel returns nothing when absent" begin
    db = make_coding_indels_db([
        ("strainA", "T1", 10, -3, "hom", "ATGC", "A"),
    ])
    @test isnothing(lookup_cds_indel(db, "T1", "strainA", 99))
    @test isnothing(lookup_cds_indel(db, "T1", "strainB", 10))
end

# ---------------------------------------------------------------------------
# apply_indel_to_cds
# ---------------------------------------------------------------------------

@testset "apply_indel_to_cds deletion" begin
    # CDS = "ATGCCCGGG", deletion at pos 4 (1-based): REF=CCCG, ALT=C
    cds = "ATGCCCGGG"
    result = apply_indel_to_cds(cds, 4, "CCCG", "C")
    @test result == "ATGCGG"
end

@testset "apply_indel_to_cds insertion" begin
    # CDS = "ATGCCC", insertion at pos 4: REF=C, ALT=CTTT
    cds = "ATGCCC"
    result = apply_indel_to_cds(cds, 4, "C", "CTTT")
    @test result == "ATGCTTTCC"
end

@testset "apply_indel_to_cds at start of sequence" begin
    cds = "ATGCCC"
    result = apply_indel_to_cds(cds, 1, "ATG", "A")
    @test result == "ACCC"
end

@testset "apply_indel_to_cds at end of sequence" begin
    cds = "ATGCCC"
    result = apply_indel_to_cds(cds, 4, "CCC", "C")
    @test result == "ATGC"
end

# ---------------------------------------------------------------------------
# Variation struct has is_processed_indel field
# ---------------------------------------------------------------------------

@testset "Variation default constructor has is_processed_indel=0" begin
    v = Variation()
    @test v.is_processed_indel == 0
end

# ---------------------------------------------------------------------------
# annotate_variations! indel path (integration)
# ---------------------------------------------------------------------------

@testset "annotate_variations! hom deletion sets product from ALT-applied CDS" begin
    # CDS = "ATGCATGAT" (9 bases)
    # Hom deletion at CDS pos 4: REF="CAT", ALT="C" → CDS becomes "ATGCGAT"
    # pic = position_in_codon(4) = 1, codon_start = 4-1 = 3 → codon = alt_seq[4:6]
    # alt_seq = "ATG" + "C" + "GAT" = "ATGCGAT" → alt_seq[4:6] = "CGA" → "R"
    cds_db = SQLite.DB()
    execute(cds_db, "CREATE TABLE coding_sequences (strain TEXT, transcript_id TEXT, sequence TEXT)")
    execute(cds_db, "INSERT INTO coding_sequences VALUES ('strainA','T1','ATGCATGAT')")
    execute(cds_db, "INSERT INTO coding_sequences VALUES ('reference','T1','ATGCATGAT')")

    indel_db = make_coding_indels_db([
        ("strainA", "T1", 4, -2, "hom", "CAT", "C"),
    ])

    tinfo = TranscriptInfo("chr1", 1, CDSInterval[])
    fs_info = precompute_frameshifts(indel_db, Dict("T1" => tinfo))

    ctx = ProcessingContext(
        "reference", Set{String}(), CDSInterval[], Dict("T1" => tinfo),
        cds_db, indel_db, fs_info, ["strainA"]
    )

    cache = TranscriptSequenceCache(Dict())
    load_transcript_if_needed!(cache, cds_db, "T1")

    annotation = PositionAnnotation(1, "T1", 1, 4, 1, "CAT", "H")

    v = Variation()
    v.strain    = "strainA"
    v.reference = "CAT"
    v.base      = "C"

    annotate_variations!([v], annotation, ctx, cache)

    @test v.is_processed_indel == 1
    @test length(v.product) == 1
    @test v.product[1] == "R"   # translate_codon("CGA") == "R"
end

@testset "annotate_variations! het insertion sets two products" begin
    # CDS = "ATGCATGAT"
    # Het insertion at CDS pos 4: REF="C", ALT="CTT"
    # REF track: extract_codon("ATGCATGAT", 4) → pic=1, codon_start=3 → "CAT" → "H"
    # ALT track: apply_indel("ATGCATGAT", 4, "C", "CTT") = "ATG" + "CTT" + "ATGAT" = "ATGCTTATGAT"
    #   pic=1, codon_start=3 → "CTT" → "L"
    cds_db = SQLite.DB()
    execute(cds_db, "CREATE TABLE coding_sequences (strain TEXT, transcript_id TEXT, sequence TEXT)")
    execute(cds_db, "INSERT INTO coding_sequences VALUES ('strainA','T1','ATGCATGAT')")
    execute(cds_db, "INSERT INTO coding_sequences VALUES ('reference','T1','ATGCATGAT')")

    indel_db = make_coding_indels_db([
        ("strainA", "T1", 4, 2, "het", "C", "CTT"),
    ])

    tinfo = TranscriptInfo("chr1", 1, CDSInterval[])
    fs_info = precompute_frameshifts(indel_db, Dict("T1" => tinfo))

    ctx = ProcessingContext(
        "reference", Set{String}(), CDSInterval[], Dict("T1" => tinfo),
        cds_db, indel_db, fs_info, ["strainA"]
    )

    cache = TranscriptSequenceCache(Dict())
    load_transcript_if_needed!(cache, cds_db, "T1")

    annotation = PositionAnnotation(1, "T1", 1, 4, 1, "CAT", "H")

    v = Variation()
    v.strain    = "strainA"
    v.reference = "C"
    v.base      = "CTT"

    annotate_variations!([v], annotation, ctx, cache)

    @test v.is_processed_indel == 1
    @test length(v.product) == 2
    @test "H" in v.product
    @test "L" in v.product   # translate_codon("CTT") == "L"
end
