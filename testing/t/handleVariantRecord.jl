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
# substitution_hgvs — Phase 1 HGVS for coding substitutions
# ---------------------------------------------------------------------------

@testset "substitution_hgvs missense: c. and p. from strand-oriented codons" begin
    (c, p) = substitution_hgvs(958, "GAC", "CAC", 1, "D", "H")
    @test c == "c.958G>C"
    @test p == "p.Asp320His"
end

@testset "substitution_hgvs synonymous uses '=' protein form" begin
    (c, p) = substitution_hgvs(402, "ACT", "ACC", 3, "T", "T")
    @test c == "c.402T>C"
    @test p == "p.Thr134="
end

@testset "substitution_hgvs nonsense uses Ter" begin
    (c, p) = substitution_hgvs(100, "CAA", "TAA", 1, "Q", "*")
    @test c == "c.100C>T"
    @test p == "p.Gln34Ter"
end

@testset "substitution_hgvs start-loss renders p.Met1?" begin
    (c, p) = substitution_hgvs(1, "ATG", "ACG", 2, "M", "T")
    @test c == "c.1T>C"
    @test p == "p.Met1?"
end

@testset "substitution_hgvs returns dots for ambiguous codon" begin
    (c, p) = substitution_hgvs(10, "ATG", "ANG", 2, "M", "X")
    @test c == "."
    @test p == "."
end

@testset "substitution_hgvs returns dot p. for unknown amino acid" begin
    (c, p) = substitution_hgvs(10, "ATG", "ACG", 2, "M", "X")
    @test c == "c.10T>C"
    @test p == "."
end

@testset "substitution_hgvs uppercases soft-masked bases" begin
    (c, p) = substitution_hgvs(958, "gac", "cac", 1, "D", "H")
    @test c == "c.958G>C"
    @test p == "p.Asp320His"
end

@testset "substitution_hgvs stop-loss emits dot p. (out of Phase 1 scope)" begin
    # ref stop(*) -> alt Gln(Q); protpos = div(979-1,3)+1 = 327
    (c, p) = substitution_hgvs(979, "TAA", "CAA", 1, "*", "Q")
    @test c == "c.979T>C"
    @test p == "."
end

# ---------------------------------------------------------------------------
# genomic_hgvs — per-allele g. (sub / del / ins / delins), left-aligned
# ---------------------------------------------------------------------------

@testset "genomic_hgvs substitution" begin
    @test genomic_hgvs("LmjF.01", 3745, "G", "C") == "LmjF.01:g.3745G>C"
end
@testset "genomic_hgvs substitution embedded in shared prefix" begin
    @test genomic_hgvs("LmjF.01", 100, "CA", "CG") == "LmjF.01:g.101A>G"
end
@testset "genomic_hgvs pure deletion, single base" begin
    @test genomic_hgvs("LmjF.01", 2531, "CA", "C") == "LmjF.01:g.2532delA"
end
@testset "genomic_hgvs pure deletion, multi base" begin
    @test genomic_hgvs("LmjF.01", 3037, "CGA", "C") == "LmjF.01:g.3038_3039delGA"
end
@testset "genomic_hgvs pure insertion" begin
    @test genomic_hgvs("LmjF.01", 1585, "C", "CGA") == "LmjF.01:g.1585_1586insGA"
end
@testset "genomic_hgvs delins from complex allele" begin
    @test genomic_hgvs("LmjF.01", 100, "ATA", "ACCTG") == "LmjF.01:g.101_102delinsCCTG"
end
@testset "genomic_hgvs delins from equal-length MNV" begin
    @test genomic_hgvs("LmjF.01", 100, "CA", "GT") == "LmjF.01:g.100_101delinsGT"
end
@testset "genomic_hgvs single-base delins" begin
    # ref core "C" (len 1), alt core "GT" -> single-position delins
    @test genomic_hgvs("LmjF.01", 100, "AC", "AGT") == "LmjF.01:g.101delinsGT"
end
@testset "genomic_hgvs uppercases soft-masked bases" begin
    @test genomic_hgvs("LmjF.01", 3745, "g", "c") == "LmjF.01:g.3745G>C"
end
@testset "genomic_hgvs returns dot for no change" begin
    @test genomic_hgvs("LmjF.01", 100, "A", "A") == "."
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
    @test ctx.sample_id_map["refStrain"] == 4
    @test length(ctx.sample_id_map) == 4
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

@testset "write_allele_file emits allele_frequency and strain_ids" begin
    sample_id_map = Dict{String,Int}("s1" => 1, "s2" => 2, "s3" => 3)

    # Three haploid hom variations: s1->A, s2->G, s3->G
    # total = 3; freq(A) = 1/3, freq(G) = 2/3
    v1 = Variation(); v1.strain = "s1"; v1.base = "A"; v1.reference = "A"; v1.coverage = "30"; v1.percent = "100"
    v2 = Variation(); v2.strain = "s2"; v2.base = "G"; v2.reference = "A"; v2.coverage = "30"; v2.percent = "100"
    v3 = Variation(); v3.strain = "s3"; v3.base = "G"; v3.reference = "A"; v3.coverage = "30"; v3.percent = "100"

    allele_buf = IOBuffer()
    write_allele_file(allele_buf, [v1, v2, v3], "LmjF.01", 100, sample_id_map)

    lines = filter(!isempty, split(String(take!(allele_buf)), "\n"))
    @test length(lines) == 2

    a_row = findfirst(l -> split(l, "\t")[3] == "A", lines)
    g_row = findfirst(l -> split(l, "\t")[3] == "G", lines)
    @test !isnothing(a_row)
    @test !isnothing(g_row)

    a_fields = split(lines[a_row], "\t")
    g_fields = split(lines[g_row], "\t")

    @test a_fields[4] == "1"         # distinct_strain_count
    @test a_fields[5] == "0.3333"    # allele_frequency: 1/3
    @test a_fields[8] == "{1}"
    @test g_fields[4] == "2"
    @test g_fields[5] == "0.6667"    # allele_frequency: 2/3
    @test g_fields[8] == "{2,3}"
end

@testset "write_allele_file allele_frequency reflects ploidy for diploid hom calls" begin
    sample_id_map = Dict{String,Int}("s1" => 1, "s2" => 2, "s3" => 3)

    # Diploid hom: s1=AA, s2=GG, s3=GG; total = 6; freq(A) = 2/6, freq(G) = 4/6
    v1 = Variation(); v1.strain = "s1"; v1.base = "A"; v1.reference = "A"; v1.coverage = "30"; v1.percent = "100"; v1.ploidy = 2
    v2 = Variation(); v2.strain = "s2"; v2.base = "G"; v2.reference = "A"; v2.coverage = "30"; v2.percent = "100"; v2.ploidy = 2
    v3 = Variation(); v3.strain = "s3"; v3.base = "G"; v3.reference = "A"; v3.coverage = "30"; v3.percent = "100"; v3.ploidy = 2

    allele_buf = IOBuffer()
    write_allele_file(allele_buf, [v1, v2, v3], "LmjF.01", 100, sample_id_map)

    lines = filter(!isempty, split(String(take!(allele_buf)), "\n"))
    a_row = findfirst(l -> split(l, "\t")[3] == "A", lines)
    g_row = findfirst(l -> split(l, "\t")[3] == "G", lines)
    @test split(lines[a_row], "\t")[5] == "0.3333"   # 2/6
    @test split(lines[g_row], "\t")[5] == "0.6667"   # 4/6
end

@testset "write_allele_file het components each have frequency 0.5" begin
    sample_id_map = Dict{String,Int}("s1" => 1)

    # Diploid het: s1 is A/G; each component weight=1, total=2, each freq=0.5
    v1 = Variation(); v1.strain = "s1"; v1.base = "A"; v1.reference = "A"; v1.alt_allele = "G"
    v1.coverage = "30"; v1.percent = "40"; v1.ploidy = 2

    allele_buf = IOBuffer()
    write_allele_file(allele_buf, [v1], "LmjF.01", 100, sample_id_map)

    lines = filter(!isempty, split(String(take!(allele_buf)), "\n"))
    a_row = findfirst(l -> split(l, "\t")[3] == "A", lines)
    g_row = findfirst(l -> split(l, "\t")[3] == "G", lines)
    @test split(lines[a_row], "\t")[5] == "0.5000"
    @test split(lines[g_row], "\t")[5] == "0.5000"
end

@testset "write_allele_file matches_reference=1 for ref allele, 0 for alt" begin
    sample_id_map = Dict{String,Int}("s1" => 1, "s2" => 2)
    v1 = Variation(); v1.strain = "s1"; v1.base = "A"; v1.reference = "A"; v1.coverage = "30"; v1.percent = "100"
    v2 = Variation(); v2.strain = "s2"; v2.base = "G"; v2.reference = "A"; v2.coverage = "30"; v2.percent = "100"

    buf = IOBuffer()
    write_allele_file(buf, [v1, v2], "LmjF.01", 100, sample_id_map)
    lines = filter(!isempty, split(String(take!(buf)), "\n"))

    a_row = findfirst(l -> split(l, "\t")[3] == "A", lines)
    g_row = findfirst(l -> split(l, "\t")[3] == "G", lines)
    @test split(lines[a_row], "\t")[9] == "1"   # ref allele
    @test split(lines[g_row], "\t")[9] == "0"   # alt allele
end

@testset "write_allele_file het: ref component matches_reference=1, alt component=0" begin
    sample_id_map = Dict{String,Int}("s1" => 1)
    v1 = Variation(); v1.strain = "s1"; v1.base = "A"; v1.reference = "A"; v1.alt_allele = "G"
    v1.coverage = "30"; v1.percent = "40"; v1.ploidy = 2

    buf = IOBuffer()
    write_allele_file(buf, [v1], "LmjF.01", 100, sample_id_map)
    lines = filter(!isempty, split(String(take!(buf)), "\n"))

    a_row = findfirst(l -> split(l, "\t")[3] == "A", lines)
    g_row = findfirst(l -> split(l, "\t")[3] == "G", lines)
    @test split(lines[a_row], "\t")[9] == "1"
    @test split(lines[g_row], "\t")[9] == "0"
end

@testset "write_allele_file has 10 columns per row" begin
    sample_id_map = Dict{String,Int}("s1" => 1)
    v1 = Variation(); v1.strain = "s1"; v1.base = "A"; v1.reference = "A"; v1.coverage = "30"; v1.percent = "100"

    buf = IOBuffer()
    write_allele_file(buf, [v1], "LmjF.01", 100, sample_id_map)
    lines = filter(!isempty, split(String(take!(buf)), "\n"))
    @test length(split(lines[1], "\t")) == 10
end

@testset "write_allele_file emits genomic_hgvs column for alt alleles" begin
    v_ref = Variation(); v_ref.strain="ref"; v_ref.base="C"; v_ref.reference="C"; v_ref.ploidy=1; v_ref.coverage="10"; v_ref.percent="100"
    v_s1  = Variation(); v_s1.strain="s1";  v_s1.base="T"; v_s1.reference="C"; v_s1.ploidy=1; v_s1.coverage="12"; v_s1.percent="100"
    v_s2  = Variation(); v_s2.strain="s2";  v_s2.base="C"; v_s2.reference="C"; v_s2.ploidy=1; v_s2.coverage="9";  v_s2.percent="100"

    buf = IOBuffer()
    write_allele_file(buf, [v_ref, v_s1, v_s2], "LmjF.01", 700, Dict("ref"=>5,"s1"=>1,"s2"=>2))
    rows = Dict{String,Vector{SubString{String}}}()
    for ln in filter(!isempty, split(String(take!(buf)), "\n"))
        f = split(ln, "\t"); rows[f[3]] = f
    end
    @test length(rows["T"]) == 10
    @test rows["T"][10] == "LmjF.01:g.700C>T"
    @test rows["C"][10] == "."
end

@testset "write_allele_file genomic_hgvs for a deletion allele" begin
    v_ref = Variation(); v_ref.strain="ref"; v_ref.base="CA"; v_ref.reference="CA"; v_ref.ploidy=1; v_ref.coverage="10"; v_ref.percent="100"
    v_s1  = Variation(); v_s1.strain="s1";  v_s1.base="C";  v_s1.reference="CA"; v_s1.ploidy=1; v_s1.coverage="12"; v_s1.percent="100"

    buf = IOBuffer()
    write_allele_file(buf, [v_ref, v_s1], "LmjF.01", 2531, Dict("ref"=>5,"s1"=>1))
    rows = Dict{String,Vector{SubString{String}}}()
    for ln in filter(!isempty, split(String(take!(buf)), "\n"))
        f = split(ln, "\t"); rows[f[3]] = f
    end
    @test rows["C"][10] == "LmjF.01:g.2532delA"
end

@testset "write_allele_file genomic_hgvs is '.' on multi-ref collision" begin
    # same allele string "C" arising from two different deletions (CA->C, CAT->C)
    # at one locus: ref span is ambiguous, so g. must be "." rather than a wrong string.
    v1 = Variation(); v1.strain="s1"; v1.base="C"; v1.reference="CA";  v1.ploidy=1; v1.coverage="10"; v1.percent="100"
    v2 = Variation(); v2.strain="s2"; v2.base="C"; v2.reference="CAT"; v2.ploidy=1; v2.coverage="11"; v2.percent="100"

    buf = IOBuffer()
    write_allele_file(buf, [v1, v2], "LmjF.01", 500, Dict("s1"=>1,"s2"=>2))
    rows = Dict{String,Vector{SubString{String}}}()
    for ln in filter(!isempty, split(String(take!(buf)), "\n"))
        f = split(ln, "\t"); rows[f[3]] = f
    end
    @test rows["C"][10] == "."
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

# ---------------------------------------------------------------------------
# build_cann_string / build_ref_cann_entry — hgvs_c / hgvs_p fields
# ---------------------------------------------------------------------------

@testset "build_cann_string appends hgvs_c and hgvs_p for a coding substitution" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=958,
                          pos_in_codon_val=1, ref_codon="GAC", ref_product="D")
    v   = make_variation(strain="s1", codon="CAC", product=["H"])
    entry = build_cann_string("G", "C", v, ann)     # SNP: ref len 1, alt len 1
    parts = split(entry, "|")
    @test length(parts) == 9
    @test parts[8] == "c.958G>C"
    @test parts[9] == "p.Asp320His"
end

@testset "build_cann_string emits dot hgvs for a pure indel" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                          pos_in_codon_val=1, ref_codon="ATG", ref_product="M")
    v   = make_variation(strain="s1", codon=".", product=String[])
    entry = build_cann_string("AT", "A", v, ann)    # deletion
    parts = split(entry, "|")
    @test length(parts) == 9
    @test parts[8] == "."
    @test parts[9] == "."
end

@testset "build_ref_cann_entry appends dot hgvs fields" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                          pos_in_codon_val=1, ref_codon="ATG", ref_product="M")
    entry = build_ref_cann_entry("r0", ann)
    parts = split(entry, "|")
    @test length(parts) == 9
    @test parts[8] == "."
    @test parts[9] == "."
end

# ---------------------------------------------------------------------------
# write_transcript_product
# ---------------------------------------------------------------------------

@testset "write_transcript_product skips DFS strains so their undefined codon does not expand to 64 rows" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                          pos_in_codon_val=2, ref_codon="ATG", ref_product="M")
    # v1: valid strain, non-DFS, codon ATT -> Ile
    v1 = make_variation(strain="s1", codon="ATT", product=["I"], downstream_of_frameshift=0)
    # v2: DFS strain, codon "." (undefined) — must not contribute any codon rows
    v2 = make_variation(strain="s2", codon=".", product=String[], downstream_of_frameshift=1)

    buf = IOBuffer()
    write_transcript_product(buf, [v1, v2], ann, "chr1", 100, Dict{String,Int}())
    lines = filter(!isempty, split(String(take!(buf)), "\n"))

    @test length(lines) == 1
    fields = split(lines[1], "\t")
    @test length(fields) == 12
    @test fields[6] == "ATT"
end

@testset "write_transcript_product deduplicates identical codons from multiple non-DFS strains" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                          pos_in_codon_val=2, ref_codon="ATG", ref_product="M")
    v1 = make_variation(strain="s1", codon="ATT", product=["I"], downstream_of_frameshift=0)
    v2 = make_variation(strain="s2", codon="ATT", product=["I"], downstream_of_frameshift=0)

    buf = IOBuffer()
    write_transcript_product(buf, [v1, v2], ann, "chr1", 100, Dict{String,Int}())
    lines = filter(!isempty, split(String(take!(buf)), "\n"))

    @test length(lines) == 1
    fields = split(lines[1], "\t")
    @test length(fields) == 12
    @test fields[6] == "ATT"
    @test fields[8] == "2"   # count = 2 strains with Ile product
end

@testset "write_transcript_product matches_ref_codon=1 when codon equals ref_codon" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                          pos_in_codon_val=1, ref_codon="ATG", ref_product="M")
    v1 = make_variation(strain="s1", codon="ATG", product=["M"], downstream_of_frameshift=0)

    buf = IOBuffer()
    write_transcript_product(buf, [v1], ann, "chr1", 100, Dict{String,Int}())
    lines = filter(!isempty, split(String(take!(buf)), "\n"))

    @test length(lines) == 1
    fields = split(lines[1], "\t")
    @test fields[10] == "1"   # matches_ref_codon
    @test fields[11] == "1"   # matches_ref_product (same product)
end

@testset "write_transcript_product matches_ref_codon=0, matches_ref_product=0 for missense codon" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                          pos_in_codon_val=2, ref_codon="ATG", ref_product="M")
    v1 = make_variation(strain="s1", codon="ATT", product=["I"], downstream_of_frameshift=0)

    buf = IOBuffer()
    write_transcript_product(buf, [v1], ann, "chr1", 100, Dict{String,Int}())
    lines = filter(!isempty, split(String(take!(buf)), "\n"))

    fields = split(lines[1], "\t")
    @test fields[10] == "0"   # matches_ref_codon: ATT != ATG
    @test fields[11] == "0"   # matches_ref_product: I != M
end

@testset "write_transcript_product matches_ref_codon=0, matches_ref_product=1 for synonymous codon" begin
    # ref_codon=TTT (Phe), alt codon=TTC (also Phe) — synonymous substitution
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                          pos_in_codon_val=3, ref_codon="TTT", ref_product="F")
    v1 = make_variation(strain="s1", codon="TTC", product=["F"], downstream_of_frameshift=0)

    buf = IOBuffer()
    write_transcript_product(buf, [v1], ann, "chr1", 100, Dict{String,Int}())
    lines = filter(!isempty, split(String(take!(buf)), "\n"))

    fields = split(lines[1], "\t")
    @test fields[10] == "0"   # matches_ref_codon: TTC != TTT
    @test fields[11] == "1"   # matches_ref_product: both encode F
end

@testset "write_transcript_product has 12 columns per row" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                          pos_in_codon_val=1, ref_codon="ATG", ref_product="M")
    v1 = make_variation(strain="s1", codon="ATT", product=["I"], downstream_of_frameshift=0)

    buf = IOBuffer()
    write_transcript_product(buf, [v1], ann, "chr1", 100, Dict{String,Int}())
    lines = filter(!isempty, split(String(take!(buf)), "\n"))
    @test length(split(lines[1], "\t")) == 12
end

@testset "write_transcript_product includes pos_in_cds and pos_in_protein" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                          pos_in_codon_val=1, ref_codon="ATG", ref_product="M")
    v1 = make_variation(strain="s1", codon="ATT", product=["I"], downstream_of_frameshift=0)

    buf = IOBuffer()
    write_transcript_product(buf, [v1], ann, "LmjF.01", 200, Dict{String,Int}())
    fields = split(filter(!isempty, split(String(take!(buf)), "\n"))[1], "\t")

    @test fields[1] == "LmjF.01"   # seq_id
    @test fields[2] == "200"        # location
    @test fields[3] == "T1"         # transcript_id
    @test fields[4] == "10"         # pos_in_cds
    @test fields[5] == "4"          # pos_in_protein: div(10-1,3)+1 = 4
end

@testset "write_transcript_product downstream_of_frameshift_strain_ids populated" begin
    sample_id_map = Dict{String,Int}("s1" => 1, "s2" => 2)
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=7,
                          pos_in_codon_val=1, ref_codon="ATG", ref_product="M")
    v1 = make_variation(strain="s1", codon="ATT", product=["I"], downstream_of_frameshift=0)
    v2 = make_variation(strain="s2", codon=".",   product=String[], downstream_of_frameshift=1)

    buf = IOBuffer()
    write_transcript_product(buf, [v1, v2], ann, "chr1", 100, sample_id_map)
    lines = filter(!isempty, split(String(take!(buf)), "\n"))

    # v2 has codon "." and is DFS — no rows for its codon, but DFS ID appears in every row
    @test length(lines) == 1
    fields = split(lines[1], "\t")
    @test fields[12] == "{2}"   # downstream_of_frameshift_strain_ids
end

@testset "write_transcript_product pos_in_protein boundary values" begin
    for (pic, expected_pip) in [(1, 1), (3, 1), (4, 2), (6, 2), (7, 3)]
        ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=pic,
                              pos_in_codon_val=((pic-1)%3)+1, ref_codon="ATG", ref_product="M")
        v1 = make_variation(strain="s1", codon="ATT", product=["I"], downstream_of_frameshift=0)
        buf = IOBuffer()
        write_transcript_product(buf, [v1], ann, "chr1", 100, Dict{String,Int}())
        fields = split(filter(!isempty, split(String(take!(buf)), "\n"))[1], "\t")
        @test fields[5] == string(expected_pip)
    end
end

# ---------------------------------------------------------------------------
# write_snp_feature — 14 genomic columns, called once per position
# ---------------------------------------------------------------------------

@testset "write_snp_feature emits 21 columns, no CDS fields" begin
    # ref strain: A; s1: T; s2: A (matches ref)
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "T"; v_s1.reference  = "A"; v_s1.ploidy = 1
    v_s2  = Variation(); v_s2.strain  = "s2";  v_s2.base  = "A"; v_s2.reference  = "A"; v_s2.ploidy = 1

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1, v_s2], 1, "LmjF.01", 500, "ref", ["s1", "s2"])
    lines = filter(!isempty, split(String(take!(buf)), "\n"))

    @test length(lines) == 1
    fields = split(lines[1], "\t")
    @test length(fields) == 21
    @test fields[1]  == "500"       # location
    @test fields[2]  == "LmjF.01"  # seq_id
    @test fields[3]  == "ref"       # reference_strain
    @test fields[4]  == "A"         # ref_allele
    @test fields[14] == "1"         # is_coding
end

@testset "write_snp_feature is_coding=0 for non-coding position" begin
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "T"; v_s1.reference  = "A"; v_s1.ploidy = 1

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1], 0, "LmjF.01", 200, "ref", ["s1"])
    fields = split(filter(!isempty, split(String(take!(buf)), "\n"))[1], "\t")
    @test fields[14] == "0"
end

# ---------------------------------------------------------------------------
# write_snp_feature — precompute columns (15-21):
#   variant_type, major_differs_from_reference, is_singleton,
#   het_strain_count, called_strain_count, no_call_strain_count, call_rate
# ---------------------------------------------------------------------------

@testset "write_snp_feature precompute columns for a SNV where major != reference" begin
    # ref=C; s1,s2,s3 hom alt T; s4 hom ref C. All 4 sequenced samples called.
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "C"; v_ref.reference = "C"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "T"; v_s1.reference  = "C"; v_s1.ploidy = 1
    v_s2  = Variation(); v_s2.strain  = "s2";  v_s2.base  = "T"; v_s2.reference  = "C"; v_s2.ploidy = 1
    v_s3  = Variation(); v_s3.strain  = "s3";  v_s3.base  = "T"; v_s3.reference  = "C"; v_s3.ploidy = 1
    v_s4  = Variation(); v_s4.strain  = "s4";  v_s4.base  = "C"; v_s4.reference  = "C"; v_s4.ploidy = 1

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1, v_s2, v_s3, v_s4], 0, "LmjF.01", 700, "ref",
                      ["s1", "s2", "s3", "s4"])
    fields = split(filter(!isempty, split(String(take!(buf)), "\n"))[1], "\t")

    @test fields[15] == "SNV"       # variant_type
    @test fields[16] == "1"         # major_differs_from_reference (major T != ref C)
    @test fields[17] == "0"         # is_singleton (minor allele C seen in 2 strains)
    @test fields[18] == "0"         # het_strain_count
    @test fields[19] == "4"         # called_strain_count (s1..s4)
    @test fields[20] == "0"         # no_call_strain_count
    @test fields[21] == "1.0000"    # call_rate
end

@testset "write_snp_feature call_rate excludes reference and reflects no-calls" begin
    # 4 sequenced samples; only s1 (ref-match) and s2 (het) have calls. s3,s4 are no-calls.
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "A"; v_s1.reference  = "A"; v_s1.ploidy = 1
    v_s2  = Variation(); v_s2.strain  = "s2";  v_s2.base  = "R"; v_s2.reference = "A"; v_s2.alt_allele = "G"; v_s2.ploidy = 2

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1, v_s2], 0, "LmjF.01", 900, "ref",
                      ["s1", "s2", "s3", "s4"])
    fields = split(filter(!isempty, split(String(take!(buf)), "\n"))[1], "\t")

    @test fields[15] == "SNV"       # variant_type (alt G is a substitution)
    @test fields[17] == "1"         # is_singleton (minor allele G in 1 strain)
    @test fields[18] == "1"         # het_strain_count (s2)
    @test fields[19] == "2"         # called_strain_count (s1, s2; ref excluded)
    @test fields[20] == "2"         # no_call_strain_count (s3, s4)
    @test fields[21] == "0.5000"    # call_rate (2/4)
end

@testset "write_snp_feature variant_type INDEL when only indel alleles present" begin
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "CA"; v_ref.reference = "CA"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "C";  v_s1.reference  = "CA"; v_s1.ploidy = 1
    v_s2  = Variation(); v_s2.strain  = "s2";  v_s2.base  = "CA"; v_s2.reference  = "CA"; v_s2.ploidy = 1

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1, v_s2], 0, "LmjF.01", 1000, "ref", ["s1", "s2"])
    fields = split(filter(!isempty, split(String(take!(buf)), "\n"))[1], "\t")

    @test fields[15] == "INDEL"
end

@testset "write_snp_feature variant_type MIXED when both snp and indel alleles present" begin
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "C";   v_ref.reference = "C"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "G";   v_s1.reference  = "C"; v_s1.ploidy = 1  # SNP
    v_s2  = Variation(); v_s2.strain  = "s2";  v_s2.base  = "CGA"; v_s2.reference  = "C"; v_s2.ploidy = 1  # insertion

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1, v_s2], 0, "LmjF.01", 1100, "ref", ["s1", "s2"])
    fields = split(filter(!isempty, split(String(take!(buf)), "\n"))[1], "\t")

    @test fields[15] == "MIXED"
end

# ---------------------------------------------------------------------------
# HSSS binary strain files
# ---------------------------------------------------------------------------

@testset "open_hsss_writers creates directories and strainIdToName.dat" begin
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1", "s2"], base)
    close_hsss_writers(state)

    for cutoff in [20, 40, 60, 80]
        dir = joinpath(base, "hsss_readFreq$(cutoff)")
        @test isdir(dir)
        lines = readlines(joinpath(dir, "strainIdToName.dat"))
        @test lines[1] == "1\tref"
        @test lines[2] == "2\ts1"
        @test lines[3] == "3\ts2"
    end
end

@testset "write_hsss_position writes binary record and contigIdToSourceId" begin
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1"], base)

    # ref: A, s1: T (passes all cutoffs)
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"
    v_ref.percent = "100"; v_ref.matches_reference = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "T"; v_s1.reference  = "A"
    v_s1.percent = "75"; v_s1.matches_reference = 0

    write_hsss_position!(state, [v_ref, v_s1], "ref", "LmjF.01", 200, ["s1"], Int8(77))
    close_hsss_writers(state)

    # Verify contig map written to all 4 dirs
    for cutoff in [20, 40, 60, 80]
        dir = joinpath(base, "hsss_readFreq$(cutoff)")
        contig_lines = readlines(joinpath(dir, "contigIdToSourceId.dat"))
        @test contig_lines[1] == "1\tLmjF.01"
    end

    # Verify binary record in readFreq20/referenceGenome.dat
    ref_bytes = read(joinpath(base, "hsss_readFreq20", "referenceGenome.dat"))
    @test length(ref_bytes) == 8
    @test ltoh(reinterpret(Int16, ref_bytes[1:2])[1]) == Int16(1)   # seq_index
    @test ltoh(reinterpret(Int32, ref_bytes[3:6])[1]) == Int32(200) # location
    @test reinterpret(Int8, ref_bytes[7:7])[1] == Int8(1)           # A = 1
    @test reinterpret(Int8, ref_bytes[8:8])[1] == Int8(77)          # product_code

    # Verify binary record for s1 (index 2)
    s1_bytes = read(joinpath(base, "hsss_readFreq20", "2"))
    @test length(s1_bytes) == 8
    @test reinterpret(Int8, s1_bytes[7:7])[1] == Int8(4)   # T = 4
    @test reinterpret(Int8, s1_bytes[8:8])[1] == Int8(77)
end

@testset "write_hsss_position skips position if only ref-matching strains at cutoff" begin
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1"], base)

    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"
    v_ref.percent = "100"; v_ref.matches_reference = 1
    # s1 matches reference
    v_s1 = Variation(); v_s1.strain = "s1"; v_s1.base = "A"; v_s1.reference = "A"
    v_s1.percent = "100"; v_s1.matches_reference = 1

    write_hsss_position!(state, [v_ref, v_s1], "ref", "chr1", 100, ["s1"], Int8(77))
    close_hsss_writers(state)

    # No records written anywhere
    @test filesize(joinpath(base, "hsss_readFreq20", "referenceGenome.dat")) == 0
    @test filesize(joinpath(base, "hsss_readFreq20", "2")) == 0
end

@testset "write_hsss_position writes unknown record for absent strain" begin
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1", "s2"], base)

    # Only ref and s1 present; s2 is absent
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"
    v_ref.percent = "100"; v_ref.matches_reference = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "T"; v_s1.reference  = "A"
    v_s1.percent = "60"; v_s1.matches_reference = 0

    write_hsss_position!(state, [v_ref, v_s1], "ref", "chr1", 50, ["s1", "s2"], Int8(0))
    close_hsss_writers(state)

    # s2 (index 3) gets unknown record
    s2_bytes = read(joinpath(base, "hsss_readFreq20", "3"))
    @test length(s2_bytes) == 8
    @test reinterpret(Int8, s2_bytes[7:7])[1] == Int8(0)  # unknown allele
    @test reinterpret(Int8, s2_bytes[8:8])[1] == Int8(0)  # unknown product
end

@testset "write_hsss_position filters by cutoff: s1 at 35% skips readFreq40" begin
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1"], base)

    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"
    v_ref.percent = "100"; v_ref.matches_reference = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "T"; v_s1.reference  = "A"
    v_s1.percent = "35"; v_s1.matches_reference = 0   # passes 20, fails 40

    write_hsss_position!(state, [v_ref, v_s1], "ref", "chr1", 100, ["s1"], Int8(0))
    close_hsss_writers(state)

    # readFreq20: s1 passes → record written
    @test filesize(joinpath(base, "hsss_readFreq20", "2")) == 8
    # readFreq40: s1 fails → position skipped (no record, not even unknown)
    @test filesize(joinpath(base, "hsss_readFreq40", "2")) == 0
end

@testset "write_hsss_position writes two records for het strain (one per allele)" begin
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1"], base)

    v_ref  = Variation(); v_ref.strain  = "ref"; v_ref.base  = "A"; v_ref.reference = "A"
    v_ref.percent = "100"; v_ref.matches_reference = 1
    # s1 is het: has both A and T alleles above threshold
    v_s1a  = Variation(); v_s1a.strain = "s1"; v_s1a.base = "A"; v_s1a.reference = "A"
    v_s1a.percent = "50"; v_s1a.matches_reference = 1
    v_s1t  = Variation(); v_s1t.strain = "s1"; v_s1t.base = "T"; v_s1t.reference = "A"
    v_s1t.percent = "50"; v_s1t.matches_reference = 0

    write_hsss_position!(state, [v_ref, v_s1a, v_s1t], "ref", "chr1", 100, ["s1"], Int8(0))
    close_hsss_writers(state)

    s1_bytes = read(joinpath(base, "hsss_readFreq20", "2"))
    # Two records: one for A allele (ref-matching, included because het), one for T allele
    @test length(s1_bytes) == 16
    allele1 = reinterpret(Int8, s1_bytes[7:7])[1]
    allele2 = reinterpret(Int8, s1_bytes[15:15])[1]
    @test Set([allele1, allele2]) == Set([Int8(1), Int8(4)])  # A=1, T=4
end

# ---------------------------------------------------------------------------
# Cache reader compatibility with transcript_product.dat format
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# compute_allele_weight_map — shared helper
# ---------------------------------------------------------------------------

@testset "compute_allele_weight_map returns correct weights for haploid hom variations" begin
    v1 = Variation(); v1.base = "A"; v1.reference = "A"; v1.ploidy = 1; v1.alt_allele = ""
    v2 = Variation(); v2.base = "T"; v2.reference = "A"; v2.ploidy = 1; v2.alt_allele = ""
    v3 = Variation(); v3.base = "T"; v3.reference = "A"; v3.ploidy = 1; v3.alt_allele = ""
    (weights, total) = compute_allele_weight_map([v1, v2, v3])
    @test weights["A"] == 1
    @test weights["T"] == 2
    @test total == 3
end

@testset "compute_allele_weight_map respects ploidy for diploid hom variations" begin
    v1 = Variation(); v1.base = "A"; v1.reference = "A"; v1.ploidy = 2; v1.alt_allele = ""
    v2 = Variation(); v2.base = "T"; v2.reference = "A"; v2.ploidy = 2; v2.alt_allele = ""
    (weights, total) = compute_allele_weight_map([v1, v2])
    @test weights["A"] == 2
    @test weights["T"] == 2
    @test total == 4
end

@testset "compute_allele_weight_map splits het call: ref and alt each get weight 1" begin
    v1 = Variation(); v1.base = "A"; v1.reference = "A"; v1.alt_allele = "G"; v1.ploidy = 2
    (weights, total) = compute_allele_weight_map([v1])
    @test weights["A"] == 1
    @test weights["G"] == 1
    @test total == 2
end

# ---------------------------------------------------------------------------
# fill_missing_coverage_gt — ploidy-aware GT string
# ---------------------------------------------------------------------------

@testset "fill_missing_coverage_gt fills GT=0/0 for diploid covered sample" begin
    record = make_vcf_record(pos=100, format_keys=["GT","DP"], sample_data=["./.:0"])
    cov    = make_coverage("s1", 99, 200, 42.0)
    result = fill_missing_coverage_gt(record, ["s1"], cov, 2)
    fmt = Dict(zip(record.format_keys, split(result[1], ":")))
    @test fmt["GT"] == "0/0"
end

@testset "fill_missing_coverage_gt fills GT=0 for haploid covered sample (ploidy=1)" begin
    record = make_vcf_record(pos=100, format_keys=["GT","DP"], sample_data=["./.:0"])
    cov    = make_coverage("s1", 99, 200, 42.0)
    result = fill_missing_coverage_gt(record, ["s1"], cov, 1)
    fmt = Dict(zip(record.format_keys, split(result[1], ":")))
    @test fmt["GT"] == "0"
end

# ---------------------------------------------------------------------------
# build_variations_from_record — coverage variation ploidy
# ---------------------------------------------------------------------------

@testset "build_variations_from_record coverage variation uses ploidy=2 for diploid" begin
    record = make_vcf_record(pos=100, format_keys=["GT","DP"], sample_data=["./.:0"])
    cov    = make_coverage("s1", 99, 200, 30.0)
    vars   = build_variations_from_record(record, ["s1"], Set{String}(), cov, 2)
    @test length(vars) == 1
    @test vars[1].ploidy == 2
    @test vars[1].matches_reference == 1
end

@testset "build_variations_from_record coverage variation defaults to ploidy=1" begin
    record = make_vcf_record(pos=100, format_keys=["GT","DP"], sample_data=["./.:0"])
    cov    = make_coverage("s1", 99, 200, 30.0)
    vars   = build_variations_from_record(record, ["s1"], Set{String}(), cov)
    @test length(vars) == 1
    @test vars[1].ploidy == 1
end

# ---------------------------------------------------------------------------
# initialize_processing_context — reads ploidy from args
# ---------------------------------------------------------------------------

@testset "initialize_processing_context reads ploidy from args" begin
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
        "ploidy"             => "2",
    )
    ctx = initialize_processing_context(args, ["strainA"])
    @test ctx.ploidy == 2
end

@testset "initialize_processing_context defaults ploidy to 1 when not in args" begin
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
    ctx = initialize_processing_context(args, ["strainA"])
    @test ctx.ploidy == 1
end

# ---------------------------------------------------------------------------
@testset "parse_cache_tsv_record reads first 4 cols of transcript_product.dat line" begin
    # transcript_product.dat has 12 tab-separated columns:
    # seq_id, location, transcript_id, pos_in_cds, pos_in_protein, codon,
    # pos_in_codon, count, product, matches_ref_codon, matches_ref_product,
    # downstream_of_frameshift_strain_ids
    line = "LmjF.01\t1234\tLmjF.01.0010\t42\t14\tATG\t1\t2\tM\t1\t1\t"
    rec = parse_cache_tsv_record(line)
    @test rec[1] == "LmjF.01"
    @test rec[2] == 1234
    @test rec[3] == "LmjF.01.0010"
    @test rec[4] == 42
end
