#!/usr/bin/env julia

# Tests for the helper functions extracted from handle_variant_record!
# Run with: julia testing/t/handleVariantRecord.jl

using Test

include(joinpath(@__DIR__, "../../bin/processSequenceVariations.jl"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function make_vcf_record(; chrom="chr1", pos=100, ref="A", alts=["T"],
                          format_keys=["GT","DP"], sample_data=["1/1:30"])
    VCFRecord(chrom, pos, ref, alts, ".", format_keys, sample_data)
end

function make_annotation(; is_coding=1, transcript_id="T1", cds_number=1,
                          pos_in_cds=42, pos_in_codon_val=2,
                          ref_codon="ATG", ref_product="M")
    PositionAnnotation(is_coding, transcript_id, cds_number, pos_in_cds,
                       pos_in_codon_val, ref_codon, ref_product)
end

# ---------------------------------------------------------------------------
# pick_snp_record
# ---------------------------------------------------------------------------

@testset "pick_snp_record returns SNP record when mixed with indels" begin
    snp   = make_vcf_record(ref="A",  alts=["T"])        # SNP: ref=1, alt=1
    indel = make_vcf_record(ref="AT", alts=["A"])         # deletion
    result = pick_snp_record([indel, snp])
    @test result.ref == "A"
    @test result.alts == ["T"]
end

@testset "pick_snp_record returns first record when only indels" begin
    del1 = make_vcf_record(ref="AT",  alts=["A"])
    del2 = make_vcf_record(ref="ATG", alts=["A"])
    result = pick_snp_record([del1, del2])
    @test result.ref == "AT"
end

@testset "pick_snp_record returns single SNP record" begin
    snp = make_vcf_record(ref="C", alts=["G"])
    @test pick_snp_record([snp]).ref == "C"
end

# ---------------------------------------------------------------------------
# build_ref_cann_entries
# ---------------------------------------------------------------------------

@testset "build_ref_cann_entries skips non-coding annotations" begin
    ann = make_annotation(is_coding=0)
    (keys, entries) = build_ref_cann_entries([ann])
    @test isempty(keys)
    @test isempty(entries)
end

@testset "build_ref_cann_entries builds r0 entry for one coding annotation" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                          pos_in_codon_val=1, ref_codon="ATG", ref_product="M")
    (keys, entries) = build_ref_cann_entries([ann])
    @test keys   == ["r0"]
    @test length(entries) == 1
    @test startswith(entries[1], "r0|ATG|M|reference|T1|10|1")
end

@testset "build_ref_cann_entries assigns r0 r1 for two coding annotations" begin
    a1 = make_annotation(transcript_id="T1")
    a2 = make_annotation(transcript_id="T2")
    (keys, entries) = build_ref_cann_entries([a1, a2])
    @test keys == ["r0", "r1"]
    @test length(entries) == 2
end

@testset "build_ref_cann_entries skips annotation where build_ref_cann_entry returns dot" begin
    # is_coding=1 but empty codon/product → build_ref_cann_entry still builds an entry
    # is_coding=0 forces "." return from build_ref_cann_entry
    ann_nc = make_annotation(is_coding=0)
    ann_c  = make_annotation(is_coding=1, transcript_id="T2")
    (keys, entries) = build_ref_cann_entries([ann_nc, ann_c])
    @test keys == ["r0"]
    @test contains(entries[1], "T2")
end

# ---------------------------------------------------------------------------
# fill_missing_coverage_gt
# ---------------------------------------------------------------------------

function make_coverage(strain, start_, end_, dp)
    # chrom_coverage: strain -> [(start, end, dp)]  (0-based half-open intervals)
    Dict{String, Vector{Tuple{Int,Int,Float64}}}(strain => [(start_, end_, Float64(dp))])
end

@testset "fill_missing_coverage_gt leaves already-set GT unchanged" begin
    record = make_vcf_record(pos=100, format_keys=["GT","DP"], sample_data=["1/1:30"])
    cov    = make_coverage("s1", 99, 200, 50.0)
    result = fill_missing_coverage_gt(record, ["s1"], cov)
    @test result[1] == "1/1:30"
end

@testset "fill_missing_coverage_gt fills GT=0 when sample missing GT but covered" begin
    record = make_vcf_record(pos=100, format_keys=["GT","DP"], sample_data=["./.:0"])
    cov    = make_coverage("s1", 99, 200, 42.0)
    result = fill_missing_coverage_gt(record, ["s1"], cov)
    fmt = Dict(zip(record.format_keys, split(result[1], ":")))
    @test fmt["GT"] == "0"
    @test fmt["DP"] == "42"
end

@testset "fill_missing_coverage_gt leaves GT missing when sample not covered" begin
    record = make_vcf_record(pos=100, format_keys=["GT","DP"], sample_data=["./.:0"])
    cov    = Dict{String, Vector{Tuple{Int,Int,Float64}}}()  # no coverage for this strain
    result = fill_missing_coverage_gt(record, ["s1"], cov)
    @test result[1] == "./.:0"
end

@testset "fill_missing_coverage_gt handles dot GT (single dot)" begin
    record = make_vcf_record(pos=100, format_keys=["GT","DP"], sample_data=[".:0"])
    cov    = make_coverage("s1", 99, 200, 10.0)
    result = fill_missing_coverage_gt(record, ["s1"], cov)
    fmt = Dict(zip(record.format_keys, split(result[1], ":")))
    @test fmt["GT"] == "0"
end

@testset "fill_missing_coverage_gt handles multiple samples" begin
    record = make_vcf_record(pos=100, format_keys=["GT","DP"],
                              sample_data=["1/1:30", "./.:0", "./.:0"])
    cov = Dict{String, Vector{Tuple{Int,Int,Float64}}}(
        "s2" => [(99, 200, 25.0)],   # s2 covered
        # s3 not covered
    )
    result = fill_missing_coverage_gt(record, ["s1", "s2", "s3"], cov)
    @test result[1] == "1/1:30"     # s1 unchanged (already set)
    fmt2 = Dict(zip(record.format_keys, split(result[2], ":")))
    @test fmt2["GT"] == "0"          # s2 filled
    @test result[3] == "./.:0"       # s3 unchanged (not covered)
end

# ---------------------------------------------------------------------------
# assign_cann_keys
# ---------------------------------------------------------------------------

@testset "assign_cann_keys assigns k0 for single entry" begin
    # alt_strain_entries: alt -> strain -> [cann_string]
    entries = Dict("T" => Dict("s1" => ["k0|ATG|M|synonymous|T1|10|1"]))
    (alt_cann, alt_ca) = assign_cann_keys(entries, ["s1"])
    @test length(alt_cann["T"]) == 1
    @test startswith(alt_cann["T"][1], "k0|")
    @test alt_ca["T"]["s1"] == ["k0"]
end

@testset "assign_cann_keys deduplicates identical entries across strains" begin
    same_entry = "k0|ATG|M|synonymous|T1|10|1"
    entries = Dict("T" => Dict(
        "s1" => [same_entry],
        "s2" => [same_entry],
    ))
    (alt_cann, alt_ca) = assign_cann_keys(entries, ["s1", "s2"])
    # Only one unique entry → only k0 in ordered list
    @test length(alt_cann["T"]) == 1
    @test alt_ca["T"]["s1"] == ["k0"]
    @test alt_ca["T"]["s2"] == ["k0"]
end

@testset "assign_cann_keys assigns k0 k1 for two distinct entries in strain order" begin
    entries = Dict("T" => Dict(
        "s1" => ["k0|ATG|M|synonymous|T1|10|1"],
        "s2" => ["k0|ATG|L|missense|T1|10|1"],
    ))
    (alt_cann, alt_ca) = assign_cann_keys(entries, ["s1", "s2"])
    # s1 processed first → k0; s2 entry is different → k1
    @test length(alt_cann["T"]) == 2
    @test startswith(alt_cann["T"][1], "k0|")
    @test startswith(alt_cann["T"][2], "k1|")
    @test alt_ca["T"]["s1"] == ["k0"]
    @test alt_ca["T"]["s2"] == ["k1"]
end

@testset "assign_cann_keys skips dot entries" begin
    entries = Dict("T" => Dict("s1" => ["."]))
    (alt_cann, alt_ca) = assign_cann_keys(entries, ["s1"])
    @test isempty(alt_cann["T"])
    @test alt_ca["T"]["s1"] == ["."]
end

@testset "assign_cann_keys handles multiple alts independently" begin
    entries = Dict(
        "T" => Dict("s1" => ["k0|ATG|M|synonymous|T1|10|1"]),
        "G" => Dict("s1" => ["k0|ATG|L|missense|T1|10|1"]),
    )
    (alt_cann, _) = assign_cann_keys(entries, ["s1"])
    @test haskey(alt_cann, "T")
    @test haskey(alt_cann, "G")
    @test startswith(alt_cann["T"][1], "k0|")
    @test startswith(alt_cann["G"][1], "k0|")
end

# ---------------------------------------------------------------------------
# collect_cann_entries_for_annotation
# ---------------------------------------------------------------------------

function make_variation(; strain="s1", base="T", matches_reference=0,
                         transcript="T1", position_in_cds=10,
                         downstream_of_frameshift=0, codon="ATT",
                         product=["I"], is_coding=1)
    v = Variation()
    v.strain = strain
    v.base   = base
    v.matches_reference = matches_reference
    v.transcript = transcript
    v.position_in_cds = position_in_cds
    v.downstream_of_frameshift = downstream_of_frameshift
    v.codon = codon
    v.product = product
    v.is_coding = is_coding
    v
end

@testset "collect_cann_entries_for_annotation skips reference strain" begin
    ann   = make_annotation(is_coding=1)
    v_ref = make_variation(strain="ref_strain")
    record = make_vcf_record(alts=["T"], sample_data=["1/1:30"])
    result = collect_cann_entries_for_annotation(
        [v_ref], ann, record, "ref_strain", ["ref_strain"],
        Dict("ref_strain" => 1)
    )
    @test isempty(result)
end

@testset "collect_cann_entries_for_annotation skips matches_reference variations" begin
    ann = make_annotation(is_coding=1)
    v   = make_variation(strain="s1", matches_reference=1)
    record = make_vcf_record(alts=["T"], sample_data=["1/1:30"])
    result = collect_cann_entries_for_annotation(
        [v], ann, record, "ref", ["s1"], Dict("s1" => 1)
    )
    @test isempty(result)
end

@testset "initialize_processing_context builds sample_id_map from all_strains" begin
    tmp_dir = mktempdir()
    transcript_db_path = joinpath(tmp_dir, "transcripts.db")
    indel_db_path      = joinpath(tmp_dir, "indels.db")
    SQLite.DB(transcript_db_path) |> close
    let db = SQLite.DB(indel_db_path)
        SQLite.execute(db, "CREATE TABLE indels (strain TEXT, transcript_id TEXT, position INTEGER, shift_amount INTEGER)")
        close(db)
    end
    gtf_path = joinpath(tmp_dir, "empty.gtf")
    write(gtf_path, "")

    args = Dict(
        "transcript_db"      => transcript_db_path,
        "indel_db"           => indel_db_path,
        "gtf_file"           => gtf_path,
        "reference_strain"   => "refStrain",
        "undone_strains_file"=> "",
    )
    all_strains = ["strainA", "strainB", "strainC"]
    ctx = initialize_processing_context(args, all_strains)

    @test ctx.sample_id_map["strainA"] == 1
    @test ctx.sample_id_map["strainB"] == 2
    @test ctx.sample_id_map["strainC"] == 3
    @test length(ctx.sample_id_map) == 3
end

@testset "write_sample_dat writes strain_id/sample_name TSV" begin
    tmp = tempname()
    write_sample_dat(["alpha", "beta", "gamma"], tmp)
    lines = readlines(tmp)
    @test lines[1] == "strain_id\tsample_name"
    @test lines[2] == "1\talpha"
    @test lines[3] == "2\tbeta"
    @test lines[4] == "3\tgamma"
    @test length(lines) == 4
end

@testset "write_allele_and_product_files emits strain_ids array in allele.dat" begin
    sample_id_map = Dict{String,Int}("s1" => 1, "s2" => 2, "s3" => 3)

    # Three hom variations: s1→A, s2→G, s3→G
    v1 = Variation(); v1.strain = "s1"; v1.base = "A"; v1.reference = "A"; v1.coverage = "30"; v1.percent = "100"
    v2 = Variation(); v2.strain = "s2"; v2.base = "G"; v2.reference = "A"; v2.coverage = "30"; v2.percent = "100"
    v3 = Variation(); v3.strain = "s3"; v3.base = "G"; v3.reference = "A"; v3.coverage = "30"; v3.percent = "100"
    variations = [v1, v2, v3]

    annotation = PositionAnnotation(0, "", 0, 0, 0, "", "")

    allele_buf  = IOBuffer()
    product_buf = IOBuffer()
    write_allele_and_product_files(allele_buf, product_buf, variations, annotation, "LmjF.01", 100, sample_id_map)

    lines = filter(!isempty, split(String(take!(allele_buf)), "\n"))
    @test length(lines) == 2   # one row per distinct allele

    a_row = findfirst(l -> split(l, "\t")[3] == "A", lines)
    g_row = findfirst(l -> split(l, "\t")[3] == "G", lines)
    @test !isnothing(a_row)
    @test !isnothing(g_row)

    a_fields = split(lines[a_row], "\t")
    g_fields = split(lines[g_row], "\t")

    @test a_fields[8] == "{1}"
    @test g_fields[8] == "{2,3}"   # sorted
end

@testset "collect_cann_entries_for_annotation returns entry keyed by alt allele" begin
    ann    = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                              pos_in_codon_val=1, ref_codon="ATG", ref_product="M")
    v      = make_variation(strain="s1", codon="ATT", product=["I"])
    record = make_vcf_record(alts=["T"], format_keys=["GT","DP"],
                              sample_data=["0/1:30"])
    result = collect_cann_entries_for_annotation(
        [v], ann, record, "ref", ["s1"], Dict("s1" => 1)
    )
    @test haskey(result, "T")
    @test haskey(result["T"], "s1")
    @test length(result["T"]["s1"]) == 1
    @test contains(result["T"]["s1"][1], "T1")
end
