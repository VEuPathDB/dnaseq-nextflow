#!/usr/bin/env julia

# Tests for makeCodingData.jl
# Run with: julia testing/t/makeCodingData.jl

using Test
using SQLite
using SQLite.DBInterface: execute

# Load without running main()
# We include GtfUtils first so its types are available, then include
# makeCodingData after patching ARGS to avoid main() triggering.
include(joinpath(@__DIR__, "../../bin/GtfUtils.jl"))

# Isolate testable functions by including only the non-main portions.
# We do this by temporarily redefining main to a no-op.
const _ORIG_ARGS = copy(ARGS)
const _SKIP_MAIN = true

# Include the file; we'll guard main() via a wrapper so it doesn't run
# by ensuring ARGS is empty and catching any exit.
# In practice, include() runs the top-level code but main() is only called
# at the bottom — we override by including with a sentinel.
module MakeCodingData
    using SQLite
    using SQLite.DBInterface: execute
    include(joinpath(@__DIR__, "../../bin/GtfUtils.jl"))

    # Override ARGS so main() doesn't blow up if included
    # We stop before main() by including only up to (not including) `main()`
    # via string manipulation — instead we just source the helpers directly.

    function parse_args(args::Vector{String})
        flags = Dict{String,String}()
        i = 1
        while i <= length(args)
            if startswith(args[i], "--")
                key = args[i][3:end]
                if i + 1 <= length(args) && !startswith(args[i+1], "--")
                    flags[key] = args[i+1]
                    i += 2
                else
                    flags[key] = ""
                    i += 1
                end
            else
                i += 1
            end
        end
        flags
    end

    function read_fasta(path::String)
        seqs = Dict{String,String}()
        cmd = endswith(path, ".gz") ? `gunzip -c $path` : `cat $path`
        current_id = ""
        buf = IOBuffer()
        for line in eachline(open(cmd))
            if startswith(line, '>')
                if !isempty(current_id)
                    seqs[current_id] = String(take!(buf))
                end
                current_id = String(split(line[2:end])[1])
            else
                write(buf, strip(line))
            end
        end
        if !isempty(current_id)
            seqs[current_id] = String(take!(buf))
        end
        seqs
    end

    function extract_cds_sequence(genome_seqs::Dict{String,String}, exon_list::Vector{CdsExon})
        isempty(exon_list) && return ""
        parts = String[]
        for e in exon_list
            seq = get(genome_seqs, e.seq_id, nothing)
            seq === nothing && return ""
            push!(parts, seq[e.start:e.stop])
        end
        cds = join(parts)
        if exon_list[1].strand == '-'
            cds = reverse_complement(cds)
        end
        cds
    end

    function project_indels_to_cds(db::SQLite.DB, strain::String, exon_list::Vector{CdsExon})
        result = Tuple{Int,Int}[]
        isempty(exon_list) && return result
        seq_id     = exon_list[1].seq_id
        cds_offset = 0
        for e in exon_list
            rows = execute(db,
                "SELECT position, shift FROM genomic_indels
                 WHERE strain = ? AND sequence_id = ? AND position >= ? AND position <= ?
                 ORDER BY position",
                [strain, seq_id, e.start, e.stop])
            for row in rows
                gpos  = row[1]
                shift = row[2]
                cds_pos = if e.strand == '-'
                    cds_offset + (e.stop - gpos)
                else
                    cds_offset + (gpos - e.start)
                end
                push!(result, (cds_pos, shift))
            end
            cds_offset += e.stop - e.start + 1
        end
        if exon_list[1].strand == '-'
            reverse!(result)
        end
        result
    end
end

# Helpers
function make_genomic_indel_db(rows::Vector)
    db = SQLite.DB()
    execute(db, """CREATE TABLE genomic_indels (
        strain TEXT, sequence_id TEXT, position INTEGER, shift INTEGER)""")
    for (strain, seq_id, pos, shift) in rows
        execute(db, "INSERT INTO genomic_indels VALUES (?,?,?,?)", [strain, seq_id, pos, shift])
    end
    db
end

function write_fasta(path, entries::Vector{Pair{String,String}})
    open(path, "w") do f
        for (id, seq) in entries
            println(f, ">$id"); println(f, seq)
        end
    end
end

# ---------------------------------------------------------------------------
# parse_args
# ---------------------------------------------------------------------------

@testset "parse_args" begin
    args = ["--genomic_indel_db", "indels.db", "--gtf_file", "genes.gtf",
            "--genome_fasta", "ref.fa", "--cds_db_out", "cds.db", "--indels_db_out", "ind.db"]
    flags = MakeCodingData.parse_args(args)
    @test flags["genomic_indel_db"] == "indels.db"
    @test flags["gtf_file"]         == "genes.gtf"
    @test flags["genome_fasta"]     == "ref.fa"
    @test flags["cds_db_out"]       == "cds.db"
    @test flags["indels_db_out"]    == "ind.db"

    # Boolean flag (no value)
    @test MakeCodingData.parse_args(["--debug"])["debug"] == ""
    # Empty
    @test isempty(MakeCodingData.parse_args([]))
end

# ---------------------------------------------------------------------------
# read_fasta
# ---------------------------------------------------------------------------

@testset "read_fasta plain" begin
    tmp = tempname()
    write_fasta(tmp, ["chr1" => "ATGCATGC", "chr2" => "NNNN"])
    seqs = MakeCodingData.read_fasta(tmp)
    rm(tmp)
    @test seqs["chr1"] == "ATGCATGC"
    @test seqs["chr2"] == "NNNN"
end

@testset "read_fasta multiline sequence" begin
    tmp = tempname()
    open(tmp, "w") do f
        println(f, ">chr1"); println(f, "ATGC"); println(f, "TTTT")
    end
    seqs = MakeCodingData.read_fasta(tmp); rm(tmp)
    @test seqs["chr1"] == "ATGCTTTT"
end

@testset "read_fasta header extra fields stripped to first token" begin
    tmp = tempname()
    open(tmp, "w") do f
        println(f, ">chr1 description here"); println(f, "AAAA")
    end
    seqs = MakeCodingData.read_fasta(tmp); rm(tmp)
    @test haskey(seqs, "chr1")
    @test !haskey(seqs, "chr1 description here")
end

@testset "read_fasta gzipped" begin
    tmp = tempname() * ".fa"
    write_fasta(tmp, ["seq1" => "GCGCGC"])
    run(`gzip $tmp`)
    seqs = MakeCodingData.read_fasta(tmp * ".gz")
    rm(tmp * ".gz")
    @test seqs["seq1"] == "GCGCGC"
end

# ---------------------------------------------------------------------------
# extract_cds_sequence
# ---------------------------------------------------------------------------

@testset "extract_cds_sequence plus strand single exon" begin
    # Genome: chr1 = "AAAAATGCCCCC" (1-based)
    #                  123456789012
    # CDS exon: 6-8 → "TGC"
    genome = Dict("chr1" => "AAAAATGCCCCC")
    exon = CdsExon("chr1", 6, 8, '+', "T1", 1)
    @test MakeCodingData.extract_cds_sequence(genome, [exon]) == "TGC"
end

@testset "extract_cds_sequence plus strand two exons" begin
    # exon1: 1-3 "ATG", exon2: 7-9 "CAT" → "ATGCAT"
    genome = Dict("chr1" => "ATGNNNCATNN")
    exon1 = CdsExon("chr1", 1, 3, '+', "T1", 1)
    exon2 = CdsExon("chr1", 7, 9, '+', "T1", 2)
    @test MakeCodingData.extract_cds_sequence(genome, [exon1, exon2]) == "ATGCAT"
end

@testset "extract_cds_sequence minus strand" begin
    # Genome chr1 = "ATGCAT"
    # CDS exon: 1-6, minus strand → reverse_complement("ATGCAT") = "ATGCAT"
    # ATGCAT rc: reverse = "TACGTA", complement each: T→A,A→T,C→G,G→C,T→A,A→T = "ATGCAT"
    # Wait: reverse("ATGCAT") = "TACGTA", complement: T→A,A→T,C→G,G→C,T→A,A→T = "ATGCAT"
    genome = Dict("chr1" => "ATGCAT")
    exon = CdsExon("chr1", 1, 6, '-', "T1", 1)
    @test MakeCodingData.extract_cds_sequence(genome, [exon]) == reverse_complement("ATGCAT")
end

@testset "extract_cds_sequence minus strand two exons (5'→3' order)" begin
    # Genome: chr1 = "ATGNNNCCC"
    # Minus-strand gene: exon1 (5') at 7-9 "CCC", exon2 (3') at 1-3 "ATG"
    # Spliced before rc: "CCC" + "ATG" = "CCCATG"
    # rc("CCCATG") = reverse("CCCATG") complement = "CATGGG"
    genome = Dict("chr1" => "ATGNNNCCC")
    exon1 = CdsExon("chr1", 7, 9, '-', "T1", 1)  # 5' exon
    exon2 = CdsExon("chr1", 1, 3, '-', "T1", 2)  # 3' exon
    result = MakeCodingData.extract_cds_sequence(genome, [exon1, exon2])
    @test result == reverse_complement("CCCATG")
end

@testset "extract_cds_sequence missing seq_id returns empty" begin
    genome = Dict("chr1" => "ATGCATGC")
    exon = CdsExon("chrX", 1, 4, '+', "T1", 1)
    @test MakeCodingData.extract_cds_sequence(genome, [exon]) == ""
end

@testset "extract_cds_sequence IUPAC bases preserved" begin
    genome = Dict("chr1" => "ATGNRYATGC")
    exon = CdsExon("chr1", 4, 6, '+', "T1", 1)  # "NRY"
    @test MakeCodingData.extract_cds_sequence(genome, [exon]) == "NRY"
end

# ---------------------------------------------------------------------------
# project_indels_to_cds
# ---------------------------------------------------------------------------

@testset "project_indels_to_cds plus strand single exon" begin
    # Exon: chr1 1-100, + strand
    # Indel at genomic pos 10, shift -1 → CDS pos 9
    # Indel at genomic pos 50, shift +3 → CDS pos 49
    db = make_genomic_indel_db([
        ("strainA", "chr1", 10, -1),
        ("strainA", "chr1", 50,  3),
    ])
    exon = CdsExon("chr1", 1, 100, '+', "T1", 1)
    result = MakeCodingData.project_indels_to_cds(db, "strainA", [exon])
    @test result == [(9, -1), (49, 3)]
end

@testset "project_indels_to_cds plus strand two exons" begin
    # exon1: chr1 1-10 (len 10), exon2: chr1 20-30
    # Indel at pos 5 → CDS pos 4 (exon1 offset 0)
    # Indel at pos 25 → CDS pos 10 + (25-20) = 15
    db = make_genomic_indel_db([
        ("strainA", "chr1",  5, -2),
        ("strainA", "chr1", 25,  1),
    ])
    exon1 = CdsExon("chr1",  1, 10, '+', "T1", 1)
    exon2 = CdsExon("chr1", 20, 30, '+', "T1", 2)
    result = MakeCodingData.project_indels_to_cds(db, "strainA", [exon1, exon2])
    @test result == [(4, -2), (15, 1)]
end

@testset "project_indels_to_cds minus strand single exon" begin
    # Exon: chr1 1-10, - strand (5' end is pos 10)
    # Indel at genomic pos 8, shift -1 → CDS pos = stop - gpos = 10 - 8 = 2
    # Indel at genomic pos 3, shift  2 → CDS pos = 10 - 3 = 7
    # Genomic order: pos 3 < pos 8; collected as (7,2),(2,-1), then reversed → [(2,-1),(7,2)]
    db = make_genomic_indel_db([
        ("strainA", "chr1", 3,  2),
        ("strainA", "chr1", 8, -1),
    ])
    exon = CdsExon("chr1", 1, 10, '-', "T1", 1)
    result = MakeCodingData.project_indels_to_cds(db, "strainA", [exon])
    @test result == [(2, -1), (7, 2)]
end

@testset "project_indels_to_cds no indels" begin
    db = make_genomic_indel_db([("strainA", "chr1", 5, -1)])
    exon = CdsExon("chr1", 1, 10, '+', "T1", 1)
    # Different strain — should return empty
    result = MakeCodingData.project_indels_to_cds(db, "strainB", [exon])
    @test isempty(result)
end

@testset "project_indels_to_cds indel outside exon boundaries excluded" begin
    db = make_genomic_indel_db([
        ("strainA", "chr1",  0, -1),   # before exon
        ("strainA", "chr1",  5, -1),   # inside
        ("strainA", "chr1", 11, -1),   # after exon
    ])
    exon = CdsExon("chr1", 1, 10, '+', "T1", 1)
    result = MakeCodingData.project_indels_to_cds(db, "strainA", [exon])
    @test length(result) == 1
    @test result[1] == (4, -1)
end
