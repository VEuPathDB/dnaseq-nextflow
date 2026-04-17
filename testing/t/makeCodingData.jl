#!/usr/bin/env julia

# Tests for makeCodingData.jl
# Run with: julia testing/t/makeCodingData.jl

using Test
using SQLite
using SQLite.DBInterface: execute

include(joinpath(@__DIR__, "../../bin/GtfUtils.jl"))
include(joinpath(@__DIR__, "../../bin/makeCodingData.jl"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function make_genomic_indel_db(rows::Vector)
    db = SQLite.DB()
    execute(db, """CREATE TABLE genomic_indels (
        strain TEXT, sequence_id TEXT, position INTEGER, shift INTEGER,
        zygosity TEXT, ref_allele TEXT, alt_allele TEXT)""")
    for (strain, seq_id, pos, shift, zygosity, ref_allele, alt_allele) in rows
        execute(db, "INSERT INTO genomic_indels VALUES (?,?,?,?,?,?,?)",
                [strain, seq_id, pos, shift, zygosity, ref_allele, alt_allele])
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
    flags = parse_args(args)
    @test flags["genomic_indel_db"] == "indels.db"
    @test flags["gtf_file"]         == "genes.gtf"
    @test flags["genome_fasta"]     == "ref.fa"
    @test flags["cds_db_out"]       == "cds.db"
    @test flags["indels_db_out"]    == "ind.db"

    # Boolean flag (no value)
    @test parse_args(["--debug"])["debug"] == ""
    # Empty
    @test isempty(parse_args(String[]))
end

# ---------------------------------------------------------------------------
# read_fasta
# ---------------------------------------------------------------------------

@testset "read_fasta plain" begin
    tmp = tempname()
    write_fasta(tmp, ["chr1" => "ATGCATGC", "chr2" => "NNNN"])
    seqs = read_fasta(tmp); rm(tmp)
    @test seqs["chr1"] == "ATGCATGC"
    @test seqs["chr2"] == "NNNN"
end

@testset "read_fasta multiline sequence" begin
    tmp = tempname()
    open(tmp, "w") do f
        println(f, ">chr1"); println(f, "ATGC"); println(f, "TTTT")
    end
    seqs = read_fasta(tmp); rm(tmp)
    @test seqs["chr1"] == "ATGCTTTT"
end

@testset "read_fasta header extra fields stripped to first token" begin
    tmp = tempname()
    open(tmp, "w") do f
        println(f, ">chr1 description here"); println(f, "AAAA")
    end
    seqs = read_fasta(tmp); rm(tmp)
    @test haskey(seqs, "chr1")
    @test !haskey(seqs, "chr1 description here")
end

@testset "read_fasta gzipped" begin
    tmp = tempname() * ".fa"
    write_fasta(tmp, ["seq1" => "GCGCGC"])
    run(`gzip $tmp`)
    seqs = read_fasta(tmp * ".gz")
    rm(tmp * ".gz")
    @test seqs["seq1"] == "GCGCGC"
end

# ---------------------------------------------------------------------------
# extract_cds_sequence
# ---------------------------------------------------------------------------

@testset "extract_cds_sequence plus strand single exon" begin
    # Genome: chr1 = "AAAAATGCCCCC" (1-based)
    # CDS exon: 6-8 → "TGC"
    genome = Dict("chr1" => "AAAAATGCCCCC")
    exon = CdsExon("chr1", 6, 8, '+', "T1", 1)
    @test extract_cds_sequence(genome, [exon]) == "TGC"
end

@testset "extract_cds_sequence plus strand two exons" begin
    # exon1: 1-3 "ATG", exon2: 7-9 "CAT" → "ATGCAT"
    genome = Dict("chr1" => "ATGNNNCATNN")
    exon1 = CdsExon("chr1", 1, 3, '+', "T1", 1)
    exon2 = CdsExon("chr1", 7, 9, '+', "T1", 2)
    @test extract_cds_sequence(genome, [exon1, exon2]) == "ATGCAT"
end

@testset "extract_cds_sequence minus strand single exon" begin
    genome = Dict("chr1" => "ATGCAT")
    exon = CdsExon("chr1", 1, 6, '-', "T1", 1)
    @test extract_cds_sequence(genome, [exon]) == reverse_complement("ATGCAT")
end

@testset "extract_cds_sequence minus strand two exons (5'→3' order)" begin
    # Minus-strand gene: exon1 (5') at 7-9 "CCC", exon2 (3') at 1-3 "ATG"
    # Each exon is individually RC'd then concatenated: RC("CCC") + RC("ATG") = "GGG" + "CAT"
    genome = Dict("chr1" => "ATGNNNCCC")
    exon1 = CdsExon("chr1", 7, 9, '-', "T1", 1)
    exon2 = CdsExon("chr1", 1, 3, '-', "T1", 2)
    @test extract_cds_sequence(genome, [exon1, exon2]) == "GGGCAT"
end

@testset "extract_cds_sequence missing seq_id returns empty" begin
    genome = Dict("chr1" => "ATGCATGC")
    exon = CdsExon("chrX", 1, 4, '+', "T1", 1)
    @test extract_cds_sequence(genome, [exon]) == ""
end

@testset "extract_cds_sequence IUPAC bases preserved" begin
    genome = Dict("chr1" => "ATGNRYATGC")
    exon = CdsExon("chr1", 4, 6, '+', "T1", 1)   # "NRY"
    @test extract_cds_sequence(genome, [exon]) == "NRY"
end

# ---------------------------------------------------------------------------
# project_indels_to_cds
# ---------------------------------------------------------------------------

@testset "project_indels_to_cds plus strand single exon" begin
    db = make_genomic_indel_db([
        ("strainA", "chr1", 10, -1, "hom", "AT",  "A"),
        ("strainA", "chr1", 50,  3, "het", "A",   "ATTT"),
    ])
    exon = CdsExon("chr1", 1, 100, '+', "T1", 1)
    result = project_indels_to_cds(db, "strainA", [exon])
    @test length(result) == 2
    @test result[1] == (9, -1, "hom", "AT", "A")
    @test result[2] == (49, 3, "het", "A", "ATTT")
end

@testset "project_indels_to_cds plus strand two exons" begin
    db = make_genomic_indel_db([
        ("strainA", "chr1",  5, -2, "hom", "ATG", "A"),
        ("strainA", "chr1", 25,  1, "het", "A",   "AT"),
    ])
    exon1 = CdsExon("chr1",  1, 10, '+', "T1", 1)
    exon2 = CdsExon("chr1", 20, 30, '+', "T1", 2)
    result = project_indels_to_cds(db, "strainA", [exon1, exon2])
    @test result[1] == (4, -2, "hom", "ATG", "A")
    @test result[2] == (15, 1, "het", "A", "AT")
end

@testset "project_indels_to_cds minus strand single exon" begin
    db = make_genomic_indel_db([
        ("strainA", "chr1", 3,  2, "het", "A",   "ATT"),
        ("strainA", "chr1", 8, -1, "hom", "AT",  "A"),
    ])
    exon = CdsExon("chr1", 1, 10, '-', "T1", 1)
    result = project_indels_to_cds(db, "strainA", [exon])
    @test result[1] == (2, -1, "hom", "AT", "A")
    @test result[2] == (7,  2, "het", "A", "ATT")
end

@testset "project_indels_to_cds no indels for strain" begin
    db = make_genomic_indel_db([("strainA", "chr1", 5, -1, "hom", "AT", "A")])
    exon = CdsExon("chr1", 1, 10, '+', "T1", 1)
    @test isempty(project_indels_to_cds(db, "strainB", [exon]))
end

@testset "project_indels_to_cds indels outside exon boundaries excluded" begin
    db = make_genomic_indel_db([
        ("strainA", "chr1",  0, -1, "hom", "AT", "A"),
        ("strainA", "chr1",  5, -1, "hom", "AT", "A"),
        ("strainA", "chr1", 11, -1, "hom", "AT", "A"),
    ])
    exon = CdsExon("chr1", 1, 10, '+', "T1", 1)
    result = project_indels_to_cds(db, "strainA", [exon])
    @test length(result) == 1
    @test result[1][1] == 4
end
