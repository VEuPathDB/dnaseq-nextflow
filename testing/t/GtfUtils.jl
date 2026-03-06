#!/usr/bin/env julia

# Tests for GtfUtils.jl
# Run with: julia testing/t/GtfUtils.jl
# GTF format matches pfal3D7.gtf: no exon_number attribute; fields are
#   seq_id, source, feature, start, stop, score, strand, frame, attributes

using Test

include(joinpath(@__DIR__, "../../bin/GtfUtils.jl"))

# Helper: write a temp GTF and return its path
function write_gtf(lines::Vector{String})
    tmp = tempname() * ".gtf"
    open(tmp, "w") do f
        for l in lines; println(f, l); end
    end
    tmp
end

# GTF line builder matching real format (no exon_number)
gtf_cds(seq_id, start, stop, strand, transcript_id) =
    "$seq_id\tPlasmoDB\tCDS\t$start\t$stop\t.\t$strand\t0\ttranscript_id \"$transcript_id\"; gene_id \"G1\"; gene_name \"G1\";"

# ---------------------------------------------------------------------------
# reverse_complement — IUPAC
# ---------------------------------------------------------------------------

@testset "reverse_complement basic" begin
    @test reverse_complement("ATGC") == "GCAT"
    @test reverse_complement("AAAA") == "TTTT"
    @test reverse_complement("A")    == "T"
    @test reverse_complement("T")    == "A"
    @test reverse_complement("G")    == "C"
    @test reverse_complement("C")    == "G"
    @test reverse_complement("")     == ""
end

@testset "reverse_complement IUPAC ambiguity codes" begin
    @test reverse_complement("R") == "Y"   # A|G  ->  C|T
    @test reverse_complement("Y") == "R"   # C|T  ->  A|G
    @test reverse_complement("S") == "S"   # G|C  ->  G|C
    @test reverse_complement("W") == "W"   # A|T  ->  A|T
    @test reverse_complement("K") == "M"   # G|T  ->  A|C
    @test reverse_complement("M") == "K"   # A|C  ->  G|T
    @test reverse_complement("B") == "V"   # C|G|T -> A|C|G
    @test reverse_complement("V") == "B"   # A|C|G -> C|G|T
    @test reverse_complement("D") == "H"   # A|G|T -> A|C|T
    @test reverse_complement("H") == "D"   # A|C|T -> A|G|T
    @test reverse_complement("N") == "N"
end

@testset "reverse_complement lowercase" begin
    @test reverse_complement("atgc") == "gcat"
    @test reverse_complement("r")    == "y"
    @test reverse_complement("n")    == "n"
    @test reverse_complement("k")    == "m"
end

@testset "reverse_complement is self-inverse" begin
    for seq in ["ATGCRYSWKMBDHVN", "atgcryswkmbdhvn", "AATTCCGG", "NNNN"]
        @test reverse_complement(reverse_complement(seq)) == seq
    end
end

# ---------------------------------------------------------------------------
# _extract_attr
# ---------------------------------------------------------------------------

@testset "_extract_attr" begin
    # Matches real GTF attribute format
    attrs = "transcript_id \"PF3D7_0200100.1\"; gene_id \"PF3D7_0200100\"; gene_name \"PF3D7_0200100\";"
    @test _extract_attr(attrs, "transcript_id") == "PF3D7_0200100.1"
    @test _extract_attr(attrs, "gene_id")        == "PF3D7_0200100"
    @test _extract_attr(attrs, "gene_name")      == "PF3D7_0200100"
    @test _extract_attr(attrs, "exon_number")    == ""   # absent in real GTF
    @test _extract_attr(attrs, "missing")        == ""
end

# ---------------------------------------------------------------------------
# parse_gtf — no exon_number (real GTF format)
# ---------------------------------------------------------------------------

@testset "parse_gtf plus-strand two-exon transcript" begin
    # Matches PF3D7_0200100.1: two CDS exons on + strand
    tmp = write_gtf([
        gtf_cds("Pf3D7_02_v3", 25232, 29035, "+", "PF3D7_0200100.1"),
        gtf_cds("Pf3D7_02_v3", 29837, 31168, "+", "PF3D7_0200100.1"),
    ])
    intervals, by_transcript = parse_gtf(tmp); rm(tmp)

    @test length(intervals) == 2
    @test haskey(by_transcript, "PF3D7_0200100.1")
    exons = by_transcript["PF3D7_0200100.1"]
    @test length(exons) == 2

    # 5'→3' on + strand: ascending genomic order
    @test exons[1].start == 25232
    @test exons[2].start == 29837

    # exon_number assigned 1, 2 via fallback
    @test exons[1].exon_number == 1
    @test exons[2].exon_number == 2
end

@testset "parse_gtf minus-strand two-exon transcript" begin
    # Matches PF3D7_0200200.1: two CDS exons on - strand
    # Genomic: 33030-33965 and 34191-34259
    # 5'→3' on minus: descending genomic → 34191-34259 is exon 1, 33030-33965 is exon 2
    tmp = write_gtf([
        gtf_cds("Pf3D7_02_v3", 33030, 33965, "-", "PF3D7_0200200.1"),
        gtf_cds("Pf3D7_02_v3", 34191, 34259, "-", "PF3D7_0200200.1"),
    ])
    intervals, by_transcript = parse_gtf(tmp); rm(tmp)

    @test haskey(by_transcript, "PF3D7_0200200.1")
    exons = by_transcript["PF3D7_0200200.1"]
    @test length(exons) == 2

    # 5'→3' on - strand: higher genomic coord first
    @test exons[1].start == 34191
    @test exons[1].stop  == 34259
    @test exons[2].start == 33030
    @test exons[2].stop  == 33965

    # exon_number assigned via fallback
    @test exons[1].exon_number == 1
    @test exons[2].exon_number == 2
end

@testset "parse_gtf intervals sorted by (seq_id, start)" begin
    tmp = write_gtf([
        gtf_cds("Pf3D7_02_v3", 29837, 31168, "+", "T1"),
        gtf_cds("Pf3D7_02_v3", 25232, 29035, "+", "T1"),
        gtf_cds("Pf3D7_02_v3", 33030, 33965, "-", "T2"),
    ])
    intervals, _ = parse_gtf(tmp); rm(tmp)

    @test intervals[1].start == 25232
    @test intervals[2].start == 29837
    @test intervals[3].start == 33030
end

@testset "parse_gtf non-CDS lines ignored" begin
    tmp = write_gtf([
        "Pf3D7_02_v3\tPlasmoDB\texon\t25232\t29035\t.\t+\t.\ttranscript_id \"T1\"; gene_id \"G1\"; gene_name \"G1\";",
        gtf_cds("Pf3D7_02_v3", 25232, 29035, "+", "T1"),
        "# a comment line",
    ])
    intervals, _ = parse_gtf(tmp); rm(tmp)

    @test length(intervals) == 1
end

# ---------------------------------------------------------------------------
# find_cds_exon
# ---------------------------------------------------------------------------

@testset "find_cds_exon" begin
    tmp = write_gtf([
        gtf_cds("Pf3D7_02_v3", 25232, 29035, "+", "T1"),
        gtf_cds("Pf3D7_02_v3", 29837, 31168, "+", "T1"),
        gtf_cds("Pf3D7_02_v3", 33030, 33965, "-", "T2"),
    ])
    intervals, _ = parse_gtf(tmp); rm(tmp)

    # Exact start of first exon
    e = find_cds_exon(intervals, "Pf3D7_02_v3", 25232)
    @test e !== nothing && e.transcript_id == "T1"

    # Mid-exon
    e = find_cds_exon(intervals, "Pf3D7_02_v3", 27000)
    @test e !== nothing && e.start == 25232

    # Exact end of first exon
    e = find_cds_exon(intervals, "Pf3D7_02_v3", 29035)
    @test e !== nothing && e.start == 25232

    # Intronic gap between exons
    e = find_cds_exon(intervals, "Pf3D7_02_v3", 29036)
    @test e === nothing

    # Second exon
    e = find_cds_exon(intervals, "Pf3D7_02_v3", 29837)
    @test e !== nothing && e.start == 29837

    # Minus-strand exon
    e = find_cds_exon(intervals, "Pf3D7_02_v3", 33500)
    @test e !== nothing && e.transcript_id == "T2"

    # Before any exon
    e = find_cds_exon(intervals, "Pf3D7_02_v3", 100)
    @test e === nothing

    # Wrong seq_id
    e = find_cds_exon(intervals, "Pf3D7_03_v3", 27000)
    @test e === nothing
end

# ---------------------------------------------------------------------------
# position_in_cds
# ---------------------------------------------------------------------------

@testset "position_in_cds plus strand" begin
    # T1: exon1 = 25232-29035 (len 3804), exon2 = 29837-31168
    tmp = write_gtf([
        gtf_cds("Pf3D7_02_v3", 25232, 29035, "+", "T1"),
        gtf_cds("Pf3D7_02_v3", 29837, 31168, "+", "T1"),
    ])
    intervals, by_transcript = parse_gtf(tmp); rm(tmp)

    exon1 = find_cds_exon(intervals, "Pf3D7_02_v3", 25232)
    @test position_in_cds(exon1, 25232, by_transcript) == 0    # first base
    @test position_in_cds(exon1, 25233, by_transcript) == 1
    @test position_in_cds(exon1, 29035, by_transcript) == 3803  # last base of exon1

    exon2 = find_cds_exon(intervals, "Pf3D7_02_v3", 29837)
    @test position_in_cds(exon2, 29837, by_transcript) == 3804  # first base of exon2
    @test position_in_cds(exon2, 29838, by_transcript) == 3805
end

@testset "position_in_cds minus strand" begin
    # T2: exon at 34191-34259 is exon1 (5'), exon at 33030-33965 is exon2
    # exon1 len = 34259 - 34191 + 1 = 69
    tmp = write_gtf([
        gtf_cds("Pf3D7_02_v3", 33030, 33965, "-", "T2"),
        gtf_cds("Pf3D7_02_v3", 34191, 34259, "-", "T2"),
    ])
    intervals, by_transcript = parse_gtf(tmp); rm(tmp)

    # Exon1 (34191-34259): pos 34259 (5'-most base) → CDS pos 0
    exon1 = find_cds_exon(intervals, "Pf3D7_02_v3", 34259)
    @test position_in_cds(exon1, 34259, by_transcript) == 0
    @test position_in_cds(exon1, 34258, by_transcript) == 1
    @test position_in_cds(exon1, 34191, by_transcript) == 68   # last base of exon1

    # Exon2 (33030-33965): first base in 5'→3' order is 33965
    exon2 = find_cds_exon(intervals, "Pf3D7_02_v3", 33965)
    @test position_in_cds(exon2, 33965, by_transcript) == 69   # exon1 len = 69
    @test position_in_cds(exon2, 33964, by_transcript) == 70
end
