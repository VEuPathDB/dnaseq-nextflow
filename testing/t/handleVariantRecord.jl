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
# extract_codon
# ---------------------------------------------------------------------------

@testset "extract_codon returns the codon containing an in-range position" begin
    seq = "ATGCGT"
    @test extract_codon(seq, 1) == "ATG"   # first position of codon 1
    @test extract_codon(seq, 3) == "ATG"   # third position of codon 1
    @test extract_codon(seq, 4) == "CGT"   # first position of codon 2
end

@testset "extract_codon returns NNN when position runs off the top of the sequence" begin
    @test extract_codon("ATG", 5) == "NNN"
end

@testset "extract_codon returns NNN for non-positive positions" begin
    # An indel-adjusted position can be <= 0 when a deletion beginning upstream
    # in genomic space extends into the CDS, so the reference position no longer
    # exists in the strain's sequence. Must not throw a BoundsError.
    @test extract_codon("ATGCGT", 0) == "NNN"
    @test extract_codon("ATGCGT", -9) == "NNN"
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

@testset "protein_hgvs missense / synonymous / start-loss / stop-loss / unknown" begin
    @test protein_hgvs(320, "D", "H") == "p.Asp320His"
    @test protein_hgvs(134, "T", "T") == "p.Thr134="
    @test protein_hgvs(1,   "M", "T") == "p.Met1?"
    @test protein_hgvs(34,  "Q", "*") == "p.Gln34Ter"
    @test protein_hgvs(327, "*", "Q") == "."
    @test protein_hgvs(10,  "M", "X") == "."
end

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

@testset "write_allele_file keys by (ref,allele), so same-string alleles from different refs stay distinct rows" begin
    # allele string "C" arising from two different deletions (CA->C, CAT->C) at one
    # locus: since aggregation is keyed by (ref, allele), each ref span gets its own
    # row with its own correct (non-ambiguous) g.HGVS, rather than collapsing into
    # a single "." row as the old allele-string-keyed aggregation did.
    v1 = Variation(); v1.strain="s1"; v1.base="C"; v1.reference="CA";  v1.ploidy=1; v1.coverage="10"; v1.percent="100"
    v2 = Variation(); v2.strain="s2"; v2.base="C"; v2.reference="CAT"; v2.ploidy=1; v2.coverage="11"; v2.percent="100"

    buf = IOBuffer()
    write_allele_file(buf, [v1, v2], "LmjF.01", 500, Dict("s1"=>1,"s2"=>2))
    rows = [split(ln, "\t") for ln in filter(!isempty, split(String(take!(buf)), "\n")) if split(ln, "\t")[3] == "C"]
    @test length(rows) == 2
    hgvs = Set(r[10] for r in rows)
    @test hgvs == Set(["LmjF.01:g.501delA", "LmjF.01:g.501_502delAT"])
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

@testset "build_cann_string reports a determinate product from an ambiguous codon" begin
    # CGN is ambiguous in its bases but not in its translation — every CGx codon
    # is Arg, so v.product collapses to ["R"]. Suppressing on the base rather than
    # the product threw this away and made CANN disagree with the HSSS files.
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=202,
                          pos_in_codon_val=1, ref_codon="GGT", ref_product="G")
    v   = make_variation(strain="s1", codon="CGN", product=["R"])
    parts = split(build_cann_string("G", "C", v, ann), "|")
    @test parts[2] == "CGN"
    @test parts[3] == "R"           # was "." before
    @test parts[4] == "missense"    # effect is now computable too
end

@testset "build_cann_string still suppresses an undetermined product" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=202,
                          pos_in_codon_val=1, ref_codon="GGT", ref_product="G")
    # Spans multiple amino acids -> no determinate product.
    v = make_variation(strain="s1", codon="NNN",
                       product=[translate_codon(c) for c in expand_codon("NNN")])
    parts = split(build_cann_string("G", "C", v, ann), "|")
    @test parts[3] == "." && parts[4] == "."

    # Unknown amino acid (no strain sequence) is also undetermined — crucially it
    # must NOT fall through, or the effect below would read "missense" for a
    # product we cannot call.
    vx = make_variation(strain="s1", codon="NNN", product=["X"])
    partsx = split(build_cann_string("G", "C", vx, ann), "|")
    @test partsx[3] == "." && partsx[4] == "."
end

@testset "build_strain_ref_cann_entry reports a determinate ambiguous codon too" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", ref_codon="GGT", ref_product="G")
    v = Variation(); v.strain="s1"; v.codon="GGN"; v.product=["G"]; v.matches_reference=1
    parts = split(build_strain_ref_cann_entry("r0", ann, v), "|")
    @test parts[2] == "GGN"
    @test parts[3] == "G"           # GGx is all Gly
    @test parts[4] == "reference"
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
    @test length(fields) == 13
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
    @test length(fields) == 13
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

@testset "write_transcript_product has 13 columns per row" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                          pos_in_codon_val=1, ref_codon="ATG", ref_product="M")
    v1 = make_variation(strain="s1", codon="ATT", product=["I"], downstream_of_frameshift=0)

    buf = IOBuffer()
    write_transcript_product(buf, [v1], ann, "chr1", 100, Dict{String,Int}())
    lines = filter(!isempty, split(String(take!(buf)), "\n"))
    @test length(split(lines[1], "\t")) == 13
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

@testset "write_transcript_product emits hgvs_p per codon row" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=958,
                          pos_in_codon_val=1, ref_codon="GAC", ref_product="D")
    v1 = make_variation(strain="s1", codon="CAC", product=["H"], downstream_of_frameshift=0)
    buf = IOBuffer()
    write_transcript_product(buf, [v1], ann, "LmjF.01", 3745, Dict{String,Int}())
    fields = split(filter(!isempty, split(String(take!(buf)), "\n"))[1], "\t")
    @test length(fields) == 13
    @test fields[13] == "p.Asp320His"
end

@testset "write_transcript_product hgvs_p is synonymous form for reference codon" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                          pos_in_codon_val=1, ref_codon="ATG", ref_product="M")
    v1 = make_variation(strain="s1", codon="ATG", product=["M"], downstream_of_frameshift=0)
    buf = IOBuffer()
    write_transcript_product(buf, [v1], ann, "chr1", 100, Dict{String,Int}())
    fields = split(filter(!isempty, split(String(take!(buf)), "\n"))[1], "\t")
    @test fields[13] == "p.Met4="
end

# ---------------------------------------------------------------------------
# write_snp_feature — 14 genomic columns, called once per position
# ---------------------------------------------------------------------------

@testset "write_snp_feature emits 31 columns, no CDS fields" begin
    # ref strain: A; s1: T; s2: A (matches ref) -> SNP alts = {T} -> snp_major = T
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "T"; v_s1.reference  = "A"; v_s1.ploidy = 1
    v_s2  = Variation(); v_s2.strain  = "s2";  v_s2.base  = "A"; v_s2.reference  = "A"; v_s2.ploidy = 1

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1, v_s2], 1, "LmjF.01", 500, "ref", ["s1", "s2"])
    lines = filter(!isempty, split(chomp(String(take!(buf))), "\n"))

    @test length(lines) == 1
    fields = split(lines[1], '\t')
    @test length(fields) == 31
    @test fields[1]  == "500"       # location
    @test fields[2]  == "LmjF.01"   # seq_id
    @test fields[3]  == "ref"       # reference_strain
    @test fields[4]  == "1"         # is_coding
    @test fields[13] == "A"                    # snp_ref_allele
    @test fields[14] == "A"                    # snp_major_allele = reference (weight 2)
    @test fields[20] == ""                     # snp_major_genomic_hgvs blank for ref slot
    @test fields[17] == "T"                    # snp_minor_allele = the alt
    @test fields[21] == "LmjF.01:g.500A>T"     # snp_minor_genomic_hgvs
end

@testset "write_snp_feature indel major_genomic_hgvs for a deletion" begin
    # ref CA at 2531; indel alts = {C} (deletion, 1 strain)
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "CA"; v_ref.reference = "CA"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "C";  v_s1.reference  = "CA"; v_s1.ploidy = 1
    v_s2  = Variation(); v_s2.strain  = "s2";  v_s2.base  = "CA"; v_s2.reference  = "CA"; v_s2.ploidy = 1

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1, v_s2], 0, "LmjF.01", 2531, "ref", ["s1", "s2"])
    fields = split(chomp(String(take!(buf))), '\t')
    @test fields[5]  == "INDEL"                    # variant_type
    @test fields[22] == "CA"                       # indel_ref_allele
    @test fields[23] == "CA"                       # indel_major_allele = reference span (weight 2)
    @test fields[29] == ""                         # indel_major_genomic_hgvs blank for ref slot
    @test fields[26] == "C"                        # indel_minor_allele = deletion
    @test fields[30] == "LmjF.01:g.2532delA"       # indel_minor_genomic_hgvs
end

@testset "write_snp_feature allele families empty when locus is monoallelic" begin
    # only the reference allele present -> no snp/indel alleles at all
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "A"; v_s1.reference  = "A"; v_s1.ploidy = 1

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1], 0, "LmjF.01", 300, "ref", ["s1"])
    fields = split(chomp(String(take!(buf))), '\t')
    @test fields[14] == ""     # snp_major_allele empty
    @test fields[23] == ""     # indel_major_allele empty
    @test fields[20] == ""     # snp_major_genomic_hgvs empty
    @test fields[29] == ""     # indel_major_genomic_hgvs empty
end

@testset "write_snp_feature is_coding=0 for non-coding position" begin
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "T"; v_s1.reference  = "A"; v_s1.ploidy = 1

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1], 0, "LmjF.01", 200, "ref", ["s1"])
    fields = split(chomp(String(take!(buf))), '\t')
    @test fields[4] == "0"
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
    fields = split(chomp(String(take!(buf))), '\t')

    @test fields[5]  == "SNV"       # variant_type
    @test fields[12] == "0"         # het_strain_count
    @test fields[7]  == "4"         # called_strain_count (s1..s4)
    @test fields[8]  == "0"         # no_call_strain_count
    @test fields[9]  == "1.0000"    # call_rate
    @test fields[14] == "T"         # snp_major_allele
    @test fields[16] == "3"         # snp_major_allele_strain_count (s1,s2,s3)
end

@testset "write_snp_feature call_rate excludes reference and reflects no-calls" begin
    # 4 sequenced samples; only s1 (ref-match) and s2 (het) have calls. s3,s4 are no-calls.
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "A"; v_s1.reference  = "A"; v_s1.ploidy = 1
    v_s2  = Variation(); v_s2.strain  = "s2";  v_s2.base  = "R"; v_s2.reference = "A"; v_s2.alt_allele = "G"; v_s2.ploidy = 2

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1, v_s2], 0, "LmjF.01", 900, "ref",
                      ["s1", "s2", "s3", "s4"])
    fields = split(chomp(String(take!(buf))), '\t')

    @test fields[5]  == "SNV"       # variant_type (alt G is a substitution)
    @test fields[12] == "1"         # het_strain_count (s2)
    @test fields[7]  == "2"         # called_strain_count (s1, s2; ref excluded)
    @test fields[8]  == "2"         # no_call_strain_count (s3, s4)
    @test fields[9]  == "0.5000"    # call_rate (2/4)
end

@testset "write_snp_feature variant_type INDEL when only indel alleles present" begin
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "CA"; v_ref.reference = "CA"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "C";  v_s1.reference  = "CA"; v_s1.ploidy = 1
    v_s2  = Variation(); v_s2.strain  = "s2";  v_s2.base  = "CA"; v_s2.reference  = "CA"; v_s2.ploidy = 1

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1, v_s2], 0, "LmjF.01", 1000, "ref", ["s1", "s2"])
    fields = split(chomp(String(take!(buf))), '\t')

    @test fields[5] == "INDEL"
end

@testset "write_snp_feature indel_frame_effect empty for non-coding indel" begin
    # Deletion CA>C (len_diff = -1) at a non-coding position (is_coding=0).
    # Frameshift/inframe classification is only meaningful inside a CDS, so the
    # frame effect must be blank when the indel is outside CDS boundaries.
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "CA"; v_ref.reference = "CA"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "C";  v_s1.reference  = "CA"; v_s1.ploidy = 1
    v_s2  = Variation(); v_s2.strain  = "s2";  v_s2.base  = "CA"; v_s2.reference  = "CA"; v_s2.ploidy = 1

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1, v_s2], 0, "LmjF.01", 1000, "ref", ["s1", "s2"])
    fields = split(chomp(String(take!(buf))), '\t')

    @test fields[5]  == "INDEL"
    @test fields[31] == ""            # indel_frame_effect blank when non-coding
end

@testset "write_snp_feature indel_frame_effect frameshift for coding indel" begin
    # Same deletion but at a coding position (is_coding=1) → frameshift.
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "CA"; v_ref.reference = "CA"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "C";  v_s1.reference  = "CA"; v_s1.ploidy = 1
    v_s2  = Variation(); v_s2.strain  = "s2";  v_s2.base  = "CA"; v_s2.reference  = "CA"; v_s2.ploidy = 1

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1, v_s2], 1, "LmjF.01", 1000, "ref", ["s1", "s2"])
    fields = split(chomp(String(take!(buf))), '\t')

    @test fields[31] == "frameshift"
end

@testset "write_snp_feature variant_type MIXED when both snp and indel alleles present" begin
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "C";   v_ref.reference = "C"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "G";   v_s1.reference  = "C"; v_s1.ploidy = 1  # SNP
    v_s2  = Variation(); v_s2.strain  = "s2";  v_s2.base  = "CGA"; v_s2.reference  = "C"; v_s2.ploidy = 1  # insertion

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1, v_s2], 0, "LmjF.01", 1100, "ref", ["s1", "s2"])
    fields = split(chomp(String(take!(buf))), '\t')

    @test fields[5] == "MIXED"
end

@testset "write_snp_feature: reference wins minor when it outranks a rarer alt" begin
    # T=3 strains, C=1 strain, ref G=2 strains (incl. synthetic ref).
    # Ranking over {ref,alts}: major=T(3), minor=ref G(2), C(1) dropped.
    vT1 = Variation(); vT1.strain="s1"; vT1.base="T"; vT1.reference="G"; vT1.ploidy=1
    vT2 = Variation(); vT2.strain="s2"; vT2.base="T"; vT2.reference="G"; vT2.ploidy=1
    vT3 = Variation(); vT3.strain="s3"; vT3.base="T"; vT3.reference="G"; vT3.ploidy=1
    vC  = Variation(); vC.strain="s4";  vC.base="C"; vC.reference="G"; vC.ploidy=1
    vR1 = Variation(); vR1.strain="s5"; vR1.base="G"; vR1.reference="G"; vR1.ploidy=1
    vRef= Variation(); vRef.strain="ref"; vRef.base="G"; vRef.reference="G"; vRef.ploidy=1

    buf = IOBuffer()
    write_snp_feature(buf, [vT1,vT2,vT3,vC,vR1,vRef], 1, "Pf3D7_01_v3", 481838, "ref",
                      ["s1","s2","s3","s4","s5"])
    f = split(chomp(String(take!(buf))), '\t')
    @test f[13] == "G"                              # snp_ref_allele
    @test f[14] == "T"                              # snp_major_allele (weight 3)
    @test f[16] == "3"                              # snp_major_allele_strain_count
    @test f[20] == "Pf3D7_01_v3:g.481838G>T"        # snp_major_genomic_hgvs
    @test f[17] == "G"                              # snp_minor_allele = REFERENCE
    @test f[19] == "2"                              # snp_minor_allele_strain_count (s5 + ref)
    @test f[21] == ""                               # snp_minor_genomic_hgvs blank for ref slot
end

@testset "write_snp_feature: reference is major at a typical ref-majority SNP" begin
    # ref C: 3 strains (incl. synthetic ref) ; alt T: 1 strain. major=ref, minor=alt.
    vR1 = Variation(); vR1.strain="s1"; vR1.base="C"; vR1.reference="C"; vR1.ploidy=1
    vR2 = Variation(); vR2.strain="s2"; vR2.base="C"; vR2.reference="C"; vR2.ploidy=1
    vT  = Variation(); vT.strain="s3";  vT.base="T"; vT.reference="C"; vT.ploidy=1
    vRef= Variation(); vRef.strain="ref"; vRef.base="C"; vRef.reference="C"; vRef.ploidy=1

    buf = IOBuffer()
    write_snp_feature(buf, [vR1,vR2,vT,vRef], 0, "chr1", 50, "ref", ["s1","s2","s3"])
    f = split(chomp(String(take!(buf))), '\t')
    @test f[14] == "C"                 # snp_major_allele = reference
    @test f[16] == "3"                 # major strain count (s1,s2,ref)
    @test f[20] == ""                  # major hgvs blank for ref slot
    @test f[17] == "T"                 # snp_minor_allele = the alt
    @test f[19] == "1"                 # minor strain count
    @test f[21] == "chr1:g.50C>T"      # minor hgvs
end

@testset "write_snp_feature: indel reference competes for indel major/minor" begin
    # ref span ACA: 2 strains (incl synthetic ref); deletion A: 1 strain.
    vR1 = Variation(); vR1.strain="s1"; vR1.base="ACA"; vR1.reference="ACA"; vR1.ploidy=1
    vD  = Variation(); vD.strain="s2";  vD.base="A";   vD.reference="ACA"; vD.ploidy=1
    vRef= Variation(); vRef.strain="ref"; vRef.base="ACA"; vRef.reference="ACA"; vRef.ploidy=1

    buf = IOBuffer()
    write_snp_feature(buf, [vR1,vD,vRef], 0, "chr1", 200, "ref", ["s1","s2"])
    f = split(chomp(String(take!(buf))), '\t')
    @test f[22] == "ACA"               # indel_ref_allele
    @test f[23] == "ACA"               # indel_major_allele = reference span
    @test f[29] == ""                  # indel major hgvs blank for ref slot
    @test f[26] == "A"                 # indel_minor_allele = deletion
    @test occursin("del", f[30])       # indel minor hgvs is a deletion
end

@testset "write_snp_feature: class with no alts stays empty (monoallelic)" begin
    vR1 = Variation(); vR1.strain="s1"; vR1.base="A"; vR1.reference="A"; vR1.ploidy=1
    vRef= Variation(); vRef.strain="ref"; vRef.base="A"; vRef.reference="A"; vRef.ploidy=1
    buf = IOBuffer()
    write_snp_feature(buf, [vR1,vRef], 0, "chr1", 10, "ref", ["s1"])
    f = split(chomp(String(take!(buf))), '\t')
    @test f[14] == ""                  # snp_major_allele empty (no snp alt)
    @test f[17] == ""                  # snp_minor_allele empty
    @test f[23] == ""                  # indel_major_allele empty
end

@testset "write_snp_feature: indel_frame_effect reflects the alt even when reference is indel-major" begin
    # ref span CA wins indel-major (weight 2), deletion C is minor (weight 1).
    # Frame effect must still report the deletion's frameshift (len C - len CA = -1).
    vR1 = Variation(); vR1.strain="s1"; vR1.base="CA"; vR1.reference="CA"; vR1.ploidy=1
    vD  = Variation(); vD.strain="s2";  vD.base="C";  vD.reference="CA"; vD.ploidy=1
    vRef= Variation(); vRef.strain="ref"; vRef.base="CA"; vRef.reference="CA"; vRef.ploidy=1
    buf = IOBuffer()
    write_snp_feature(buf, [vR1,vD,vRef], 1, "LmjF.01", 2531, "ref", ["s1","s2"])
    f = split(chomp(String(take!(buf))), '\t')
    @test f[23] == "CA"           # indel_major_allele = reference span
    @test f[31] == "frameshift"   # frame effect still from the deletion alt
end


# ---------------------------------------------------------------------------
# HSSS binary strain files
# ---------------------------------------------------------------------------

# Helper: read the 8-byte HSSS records out of a file as (seq, loc, allele, product).
function read_hsss(path)
    bytes = read(path)
    @assert length(bytes) % 8 == 0 "HSSS file $(path) is not a whole number of 8-byte records"
    [(ltoh(reinterpret(Int16, bytes[i+1:i+2])[1]),
      ltoh(reinterpret(Int32, bytes[i+3:i+6])[1]),
      reinterpret(Int8, bytes[i+7:i+7])[1],
      reinterpret(Int8, bytes[i+8:i+8])[1])
     for i in 0:8:(length(bytes)-8)]
end

hsss_alleles_of(path)  = [r[3] for r in read_hsss(path)]
hsss_products_of(path) = [r[4] for r in read_hsss(path)]

# Build a non-reference Variation with resolved GT slots.
function hsss_var(strain, slots; reference="A", codon="", pic=2, percent="100")
    v = Variation()
    v.strain = strain
    v.reference = reference
    v.base = length(unique(slots)) == 1 ? slots[1] : "N"
    v.allele_slots = slots
    v.percent = percent
    v.matches_reference = all(==(reference), slots) ? 1 : 0
    v.codon = codon
    v.position_in_codon = pic
    v.is_coding = isempty(codon) ? 0 : 1
    # Mirror annotate_variations!: product is derived from the strain's consensus
    # codon, expanded over any IUPAC ambiguity. HSSS reads this field, so the
    # fixture has to populate it the same way the pipeline does.
    v.product = isempty(codon) ? String[] : [translate_codon(c) for c in expand_codon(codon)]
    v
end

@testset "open_hsss_writers creates a file for the reference strain too" begin
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

        # The reference IS a selectable sample, so file "1" must exist. It stays
        # empty: the reference matches itself everywhere, so it has zero records,
        # and hsssFindPolymorphic folds it in via ref_count.
        @test isfile(joinpath(dir, "1"))
        @test filesize(joinpath(dir, "1")) == 0
    end
end

@testset "product_for_allele resolves the SNP slot against the strain codon" begin
    # Codon ARG with the ambiguity at position 2: allele A -> AAG (K), allele G -> AGG (R)
    @test product_for_allele("ARG", 2, "A") == "K"
    @test product_for_allele("ARG", 2, "G") == "R"
    # Unambiguous codon: substitution still applies
    @test product_for_allele("AAG", 2, "G") == "R"
    # Non-coding / no codon
    @test product_for_allele("", 2, "A") == ""
end

@testset "product_for_allele returns X when a second site in the codon is ambiguous" begin
    # RRG: resolving position 2 still leaves position 1 ambiguous -> AAG(K) vs GAG(E)
    @test product_for_allele("RRG", 2, "A") == "X"
end

@testset "product_for_allele complements the allele on the minus strand" begin
    # Codons come from the CDS (reverse-complemented for minus-strand genes);
    # alleles are plus-strand genomic. The same genomic allele must therefore
    # translate differently depending on strand.
    @test product_for_allele("ATG", 2, "A",  1) == "K"   # plus:  A     -> AAG
    @test product_for_allele("ATG", 2, "A", -1) == "M"   # minus: A->T  -> ATG
end

@testset "product_for_allele is a no-op for a strain's own homozygous allele" begin
    # THE invariant that catches strand errors. A strain's CDS codon already
    # contains that strain's allele, so substituting the same allele back in must
    # reproduce translate_codon(codon) exactly — on either strand. Violating this
    # is what manufactured spurious stop codons on minus-strand genes.
    @test product_for_allele("ACG", 2, "G", -1) == translate_codon("ACG")
    @test product_for_allele("AGG", 2, "G",  1) == translate_codon("AGG")

    # Exhaustive: every codon, every position, both strands.
    comp = Dict('A'=>'T','T'=>'A','C'=>'G','G'=>'C')
    for b1 in "ACGT", b2 in "ACGT", b3 in "ACGT", pic in 1:3, strand in (1, -1)
        codon    = string(b1, b2, b3)
        cds_base = codon[pic]
        genomic  = strand == -1 ? string(comp[cds_base]) : string(cds_base)
        @test product_for_allele(codon, pic, genomic, strand) == translate_codon(codon)
    end
end

@testset "products_for_alleles does not substitute into an unambiguous codon" begin
    # The strain's own codon is authoritative. With an unambiguous codon every
    # allele maps to translate_codon(codon) and no substitution happens — so a
    # codon whose bases disagree with the VCF genotype (consensus filtering,
    # indel shift) can no longer invent an amino acid no strain actually has.
    @test products_for_alleles(["R"], "AGG", 2, ["G"], 1)  == ["R"]
    @test products_for_alleles(["R"], "AGG", 2, ["G"], -1) == ["R"]

    # Deliberately inconsistent input: allele "A" does not match the codon's base.
    # Old behaviour substituted anyway and returned a fabricated product; now the
    # codon wins on both strands and for any position.
    for pic in 1:3, strand in (1, -1), allele in ["A","C","G","T"]
        @test products_for_alleles(["Y"], "TAT", pic, [allele], strand) == [translate_codon("TAT")]
    end

    # Specifically: TAT + allele A at position 3 used to fabricate TAA (a stop).
    @test products_for_alleles(["Y"], "TAT", 3, ["A"], 1) == ["Y"]
    @test product_for_allele("TAT", 3, "A", 1)     == "*"   # the raw substitution still would
end

@testset "products_for_alleles resolves each allele when the codon is ambiguous" begin
    # Het: codon carries an IUPAC code, so the alleles genuinely differ and
    # substitution is the only way to tell them apart.
    @test products_for_alleles(["K","R"], "ARG", 2, ["A", "G"], 1) == ["K", "R"]
    # Minus strand: genomic A/G complement to T/C -> ATG(M), ACG(T)
    @test products_for_alleles(["M","T"], "AYG", 2, ["A", "G"], -1) == ["M", "T"]
end

@testset "products_for_alleles handles missing and frameshifted codons" begin
    # No strain sequence: coding, but the amino acid is unknown -> X, which the C
    # normalizes to -1 and ignores on the strain side.
    @test products_for_alleles(["X"], "NNN", 2, ["A"], 1) == ["X"]
    # Downstream of a frameshift, and non-coding: no product at all -> "" -> byte 0.
    # Deliberately NOT "X", which would assert "coding but unknown".
    @test products_for_alleles(String[], ".", 2, ["A"], 1) == [""]
    @test products_for_alleles(String[], "", 2, ["A"], 1) == [""]
    @test hsss_product_code("")  == Int8(0)
    @test products_for_alleles(["R"], "AGG", 2, String[], 1) == String[]
end

@testset "products_for_alleles reads Variation.product, not the codon" begin
    # Variation.product is the single source of truth — the VEuPath product call,
    # computed once from the strain's consensus codon. products_for_alleles must
    # report it rather than re-deriving from the codon, so the three consumers
    # (transcript_product.dat, CANN, HSSS) cannot drift apart. Here the product
    # list deliberately disagrees with the codon: the product wins.
    @test products_for_alleles(["W"], "AGG", 2, ["G"], 1) == ["W"]
end

# ---------------------------------------------------------------------------
# CANN reference entries: per-strain when the consensus codon differs
# ---------------------------------------------------------------------------

@testset "build_strain_ref_cann_entry reports the strain's own product" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=42,
                          pos_in_codon_val=2, ref_codon="AAG", ref_product="K")
    # Strain carries the REFERENCE allele here but another variant in the same
    # codon, so its consensus codon (GAG -> E) differs from the reference (AAG -> K).
    v = Variation(); v.strain="s1"; v.codon="GAG"; v.product=["E"]; v.matches_reference=1
    f = split(build_strain_ref_cann_entry("r0", ann, v), "|")
    @test f[1] == "r0"
    @test f[2] == "GAG"          # the strain's codon, not AAG
    @test f[3] == "E"            # the strain's product, not K
    @test f[4] == "reference"
    @test f[5] == "T1"

    # An unknown product (no strain sequence) reports nothing, mirroring the
    # alt-entry rule — "X" is not a callable amino acid.
    vn = Variation(); vn.strain="s2"; vn.codon="NNN"; vn.product=["X"]; vn.matches_reference=1
    @test split(build_strain_ref_cann_entry("r0", ann, vn), "|")[3] == "."
end

@testset "assign_ref_cann_keys appends per-strain entries and dedupes them" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", ref_codon="AAG", ref_product="K")
    v1 = Variation(); v1.codon="GAG"; v1.product=["E"]
    v2 = Variation(); v2.codon="CAG"; v2.product=["Q"]
    e1 = build_strain_ref_cann_entry("r0", ann, v1)
    e2 = build_strain_ref_cann_entry("r0", ann, v2)

    # s1 and s3 share a codon so must share a key; s2 gets its own.
    refents = Dict("s1"=>[e1], "s2"=>[e2], "s3"=>[e1])
    (ref_keys, entries, s2r) = assign_ref_cann_keys([ann], refents, ["s1","s2","s3"])

    @test ref_keys == ["r0"]                 # the shared reference entry
    @test length(entries) == 3               # r0 + two distinct strain entries
    @test s2r["s1"] == ["r1"]
    @test s2r["s3"] == ["r1"]                # deduplicated onto the same key
    @test s2r["s2"] == ["r2"]
    @test any(startswith(e, "r1|GAG|E|") for e in entries)
    @test any(startswith(e, "r2|CAG|Q|") for e in entries)
end

@testset "assign_ref_cann_keys leaves reference-codon strains on the shared entry" begin
    ann = make_annotation(is_coding=1, ref_codon="AAG", ref_product="K")
    (ref_keys, entries, s2r) = assign_ref_cann_keys([ann], Dict{String,Vector{String}}(), ["s1"])
    @test ref_keys == ["r0"]
    @test length(entries) == 1
    @test isempty(s2r)
end

@testset "build_ca_values gives a strain its own r-key when it has one" begin
    rec = make_vcf_record(pos=100, format_keys=["GT","DP"], sample_data=["0:30", "0:30"])
    ca = build_ca_values(rec.format_keys, rec.sample_data, rec.alts, "T",
                         ["r0"], ["s1","s2"],
                         Dict{String,Vector{String}}(),
                         Dict("s1"=>["r1"]))
    @test ca[1] == "r1"    # s1 -> its own product
    @test ca[2] == "r0"    # s2 -> shared reference product
end

@testset "hsss_product_code encodes ascii, empty as 0" begin
    @test hsss_product_code("K") == Int8(75)
    @test hsss_product_code("*") == Int8(42)
    @test hsss_product_code("")  == Int8(0)
end

@testset "write_hsss_position gives each strain its own product byte" begin
    # THE core regression for the collapsed-product bug: two strains whose codons
    # translate differently must NOT share a product byte, or hsssFindPolymorphic
    # can never report non-syn (it detects it only by the byte differing).
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1", "s2"], base)

    s1 = hsss_var("s1", ["G"]; codon="AGG", pic=2)   # -> R
    s2 = hsss_var("s2", ["T"]; codon="ATG", pic=2)   # -> M

    write_hsss_position!(state, [s1, s2], "A", "K", "chr1", 100, ["s1", "s2"])
    close_hsss_writers(state)

    d = joinpath(base, "hsss_readFreq20")
    @test hsss_products_of(joinpath(d, "2")) == [Int8(codepoint('R'))]
    @test hsss_products_of(joinpath(d, "3")) == [Int8(codepoint('M'))]
    @test hsss_products_of(joinpath(d, "referenceGenome.dat")) == [Int8(codepoint('K'))]
end

@testset "write_hsss_position never writes X into referenceGenome.dat" begin
    # refProduct is NOT X-normalized in hsssFindPolymorphic.c (unlike the strain
    # path), so 'X'=88 there reads as a real amino acid and fires a spurious
    # non-syn. An unknown reference product must be written as 0 (non-coding).
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1"], base)
    s1 = hsss_var("s1", ["G"]; codon="AGG", pic=2)
    write_hsss_position!(state, [s1], "A", "X", "chr1", 100, ["s1"])
    close_hsss_writers(state)

    d = joinpath(base, "hsss_readFreq20")
    @test hsss_products_of(joinpath(d, "referenceGenome.dat")) == [Int8(0)]
end

@testset "write_hsss_position writes one record per allele for a diploid het" begin
    # Regression for IUPAC compression: a het must NOT collapse to a single
    # record with an ambiguity base (which maps to allele code 0 = unknown).
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1"], base)

    s1 = hsss_var("s1", ["A", "G"]; codon="ARG", pic=2)   # het A/G, ref A
    write_hsss_position!(state, [s1], "A", "K", "chr1", 100, ["s1"])
    close_hsss_writers(state)

    d = joinpath(base, "hsss_readFreq20")
    recs = read_hsss(joinpath(d, "2"))
    @test length(recs) == 2
    # Real allele codes, not 0/unknown. Ref-matching slot included, per the Perl.
    @test Set(r[3] for r in recs) == Set([Int8(1), Int8(3)])   # A=1, G=3
    # Each allele carries its own product: A->K, G->R
    byallele = Dict(r[3] => r[4] for r in recs)
    @test byallele[Int8(1)] == Int8(codepoint('K'))
    @test byallele[Int8(3)] == Int8(codepoint('R'))
end

@testset "write_hsss_position writes one record for a homozygous diploid call" begin
    # 1/1 resolves to ["G","G"] — that is one allele, not two, and must not be
    # double-counted in the C's alleleCount.
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1"], base)
    s1 = hsss_var("s1", ["G", "G"]; codon="AGG", pic=2)
    write_hsss_position!(state, [s1], "A", "K", "chr1", 100, ["s1"])
    close_hsss_writers(state)

    @test length(read_hsss(joinpath(base, "hsss_readFreq20", "2"))) == 1
end

@testset "write_hsss_position skips position when no strain differs from reference" begin
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1"], base)
    s1 = hsss_var("s1", ["A"]; codon="AAG", pic=2)   # matches reference
    write_hsss_position!(state, [s1], "A", "K", "chr1", 100, ["s1"])
    close_hsss_writers(state)

    @test filesize(joinpath(base, "hsss_readFreq20", "referenceGenome.dat")) == 0
    @test filesize(joinpath(base, "hsss_readFreq20", "2")) == 0
end

@testset "write_hsss_position writes contigIdToSourceId and unknown for absent strain" begin
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1", "s2"], base)
    s1 = hsss_var("s1", ["T"]; codon="ATG", pic=2, percent="60")
    write_hsss_position!(state, [s1], "A", "K", "LmjF.01", 50, ["s1", "s2"])
    close_hsss_writers(state)

    for cutoff in [20, 40, 60, 80]
        d = joinpath(base, "hsss_readFreq$(cutoff)")
        @test readlines(joinpath(d, "contigIdToSourceId.dat"))[1] == "1\tLmjF.01"
    end

    # s2 absent -> unknown allele and unknown product
    recs = read_hsss(joinpath(base, "hsss_readFreq20", "3"))
    @test length(recs) == 1
    @test recs[1][3] == Int8(0)
    @test recs[1][4] == Int8(0)
end

@testset "write_hsss_position filters by cutoff: 35% passes 20 but not 40" begin
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1"], base)
    s1 = hsss_var("s1", ["T"]; codon="ATG", pic=2, percent="35")
    write_hsss_position!(state, [s1], "A", "K", "chr1", 100, ["s1"])
    close_hsss_writers(state)

    @test filesize(joinpath(base, "hsss_readFreq20", "2")) == 8
    @test filesize(joinpath(base, "hsss_readFreq40", "2")) == 0
end

@testset "select_hsss_annotation picks the transcript with the longest CDS" begin
    # Alternative transcripts of one gene give the same amino acid at a position,
    # so this only bites for overlapping genes — but it must be deterministic and
    # must never fall back to "disagreement => non-coding".
    tinfo = Dict(
        "short" => TranscriptInfo("chr1", 1, [CDSInterval("chr1", 1, 30, 1, "short", 1)]),
        "long"  => TranscriptInfo("chr1", 1, [CDSInterval("chr1", 1, 60, 1, "long", 1)]),
    )
    anns = [make_annotation(transcript_id="short", ref_product="K"),
            make_annotation(transcript_id="long",  ref_product="M")]

    @test select_hsss_annotation(anns, tinfo).transcript_id == "long"
    # Order-independent
    @test select_hsss_annotation(reverse(anns), tinfo).transcript_id == "long"
end

@testset "select_hsss_annotation prefers a coding annotation over non-coding" begin
    tinfo = Dict(
        "nc" => TranscriptInfo("chr1", 1, [CDSInterval("chr1", 1, 90, 1, "nc", 1)]),
        "c"  => TranscriptInfo("chr1", 1, [CDSInterval("chr1", 1, 30, 1, "c", 1)]),
    )
    anns = [make_annotation(transcript_id="nc", is_coding=0),
            make_annotation(transcript_id="c",  is_coding=1)]
    @test select_hsss_annotation(anns, tinfo).transcript_id == "c"
end

# ---------------------------------------------------------------------------
# Cache reader compatibility with transcript_product.dat format
# ---------------------------------------------------------------------------

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
    vars   = build_variations_from_record(record, ["s1"], cov, 2)
    @test length(vars) == 1
    @test vars[1].ploidy == 2
    @test vars[1].matches_reference == 1
end

@testset "build_variations_from_record coverage variation defaults to ploidy=1" begin
    record = make_vcf_record(pos=100, format_keys=["GT","DP"], sample_data=["./.:0"])
    cov    = make_coverage("s1", 99, 200, 30.0)
    vars   = build_variations_from_record(record, ["s1"], cov)
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

# ---------------------------------------------------------------------------
# -m both: cross-record ./. must not fabricate a duplicate reference variation
# ---------------------------------------------------------------------------

@testset "build_variations_from_record: synthesize_ref=false skips covered ./." begin
    indel = make_vcf_record(ref="ATG", alts=["A"],
                            format_keys=["GT","DP"], sample_data=["./.:0"])
    chrom_cov = Dict("S1" => [(1, 1000, 30.0)])
    vars = build_variations_from_record(indel, ["S1"], chrom_cov, 1;
                                        synthesize_ref=false)
    @test isempty(vars)
end

@testset "build_variations_from_record: synthesize_ref=true keeps old behavior" begin
    indel = make_vcf_record(ref="ATG", alts=["A"],
                            format_keys=["GT","DP"], sample_data=["./.:0"])
    chrom_cov = Dict("S1" => [(1, 1000, 30.0)])
    vars = build_variations_from_record(indel, ["S1"], chrom_cov, 1;
                                        synthesize_ref=true)
    @test length(vars) == 1
    @test vars[1].matches_reference == 1
end

@testset "build_cann_string: indel downstream of frameshift → compound effect" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=42,
                          pos_in_codon_val=2)
    v = Variation()
    v.downstream_of_frameshift = 1
    # deletion ATG>A : len_diff = -2 → structurally a frameshift
    s = build_cann_string("ATG", "A", v, ann)
    @test occursin("frameshift&downstream_frameshift", s)
end

@testset "build_cann_string: indel NOT downstream of frameshift → structural only" begin
    ann = make_annotation(is_coding=1)
    v = Variation()
    v.downstream_of_frameshift = 0
    s = build_cann_string("ATG", "A", v, ann)   # frameshift
    @test occursin("|frameshift|", s)
    @test !occursin("downstream_frameshift", s)
end

# ---------------------------------------------------------------------------
# handle_variant_record! — with bcftools -m both, a locus can carry BOTH a
# SNP class-record and an indel class-record. BOTH must contribute output.vcf
# rows (each reading GT from its OWN sample columns), not just the SNP record.
# ---------------------------------------------------------------------------

# Builds a ProcessingContext whose GTF is empty, so every position is
# intergenic (non-coding) and no transcript/indel DB queries are exercised.
function make_intergenic_ctx(all_strains::Vector{String}; reference_strain="ref", ploidy=1)
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
        "transcript_db"       => transcript_db_path,
        "indel_db"            => indel_db_path,
        "gtf_file"            => gtf_path,
        "reference_strain"    => reference_strain,
        "ploidy"              => string(ploidy),
    )
    initialize_processing_context(args, all_strains)
end

@testset "handle_variant_record! emits BOTH the SNP and indel class-record rows" begin
    all_strains = ["S1", "S2"]
    ctx = make_intergenic_ctx(all_strains)

    vcf_buf = IOBuffer()
    hsss    = open_hsss_writers("ref", all_strains, mktempdir())
    writers = OutputWriters(vcf_buf, IOBuffer(), IOBuffer(), IOBuffer(), hsss)
    transcript_cache = TranscriptSequenceCache(Dict{String, Dict{String,String}}())

    # SNP record:   S1 = A>G (1/1), S2 = no-call.
    # Indel record: S1 = no-call, S2 = ATG>A (1/1) deletion.
    # Same locus, same strain column order.
    snp   = VCFRecord("chr1", 100, "A",   ["G"], ".", ["GT","DP"], ["1/1:30", "./.:0"])
    indel = VCFRecord("chr1", 100, "ATG", ["A"], ".", ["GT","DP"], ["./.:0", "1/1:30"])

    # chrom_coverage is keyed by STRAIN name (not chromosome).
    chrom_cov = Dict{String, Vector{Tuple{Int,Int,Float64}}}(
        "S1" => [(1, 1000, 30.0)],
        "S2" => [(1, 1000, 30.0)],
    )

    handle_variant_record!([snp, indel], Tuple{String,Int}[], ctx, writers,
                           transcript_cache, all_strains, chrom_cov)
    close_hsss_writers(hsss)

    rows = [split(l, '\t') for l in filter(!isempty, split(String(take!(vcf_buf)), "\n"))]
    # (REF, ALT) pairs written for this locus
    ref_alt = Set((r[4], r[5]) for r in rows)

    @test ("A", "G") in ref_alt        # SNP class-record row
    @test ("ATG", "A") in ref_alt      # indel class-record row (dropped before the fix)
    # POS/CHROM are shared across both class-records
    @test all(r -> r[1] == "chr1" && r[2] == "100", rows)
end

# ---------------------------------------------------------------------------
# handle_variant_record! — allele.dat must not double-count a strain that is
# ./. on the SNP record precisely because it carries a REAL call on the sibling
# insertion record. Such a strain must NOT get a fabricated reference variation.
# Use an INSERTION (A>AT); a co-located SNP+deletion has a separate base-collapse
# bug that is out of scope here.
# ---------------------------------------------------------------------------

@testset "handle_variant_record! does not synthesize reference for indel-carrying strain" begin
    all_strains = ["S1", "S2"]
    ctx = make_intergenic_ctx(all_strains)

    vcf_buf    = IOBuffer()
    allele_buf = IOBuffer()
    hsss       = open_hsss_writers("ref", all_strains, mktempdir())
    # OutputWriters fields: vcf_fh, snp_fh, allele_fh, transcript_product_fh, hsss
    writers = OutputWriters(vcf_buf, IOBuffer(), allele_buf, IOBuffer(), hsss)
    transcript_cache = TranscriptSequenceCache(Dict{String, Dict{String,String}}())

    # SNP record:       S1 = A>G (1/1),  S2 = no-call.
    # Insertion record: S1 = no-call,    S2 = A>AT (1/1).
    # S2 is ./. on the SNP record ONLY because it carries the insertion.
    snp = VCFRecord("chr1", 100, "A", ["G"],  ".", ["GT","DP"], ["1/1:30", "./.:0"])
    ins = VCFRecord("chr1", 100, "A", ["AT"], ".", ["GT","DP"], ["./.:0", "1/1:30"])

    chrom_cov = Dict{String, Vector{Tuple{Int,Int,Float64}}}(
        "S1" => [(1, 1000, 30.0)],
        "S2" => [(1, 1000, 30.0)],
    )

    handle_variant_record!([snp, ins], Tuple{String,Int}[], ctx, writers,
                           transcript_cache, all_strains, chrom_cov)
    close_hsss_writers(hsss)

    # allele.dat columns: location, seq_id, allele, distinct_strain_count,
    #                     frequency, avg_coverage, avg_percent, strain_ids,
    #                     matches_reference, genomic_hgvs
    arows = [split(l, '\t') for l in filter(!isempty, split(String(take!(allele_buf)), "\n"))]
    by_allele = Dict(r[3] => r for r in arows)

    # Reference allele A must be credited with the reference strain ONLY (count 1),
    # NOT with S2 (which carries the insertion). Alt alleles G and AT each carry 1.
    @test haskey(by_allele, "A")
    @test parse(Int, by_allele["A"][4]) == 1
    @test haskey(by_allele, "G")
    @test parse(Int, by_allele["G"][4]) == 1
    @test haskey(by_allele, "AT")
    @test parse(Int, by_allele["AT"][4]) == 1

    # Invariant: no strain counted under more than one allele.
    # Total distinct_strain_count across all alleles == 3 (S1→G, S2→AT, ref→A), NOT 4.
    @test sum(parse(Int, r[4]) for r in arows) == 3
end

# ---------------------------------------------------------------------------
# aggregate_locus_alleles / classify_allele — (ref,base)-keyed aggregation
# ---------------------------------------------------------------------------

function mkvar(; strain, reference, base, alt_allele="", ploidy=1,
                coverage="30", percent="100", matches_reference=0)
    v = Variation()
    v.strain = strain; v.reference = reference; v.base = base
    v.alt_allele = alt_allele; v.ploidy = ploidy
    v.coverage = coverage; v.percent = percent
    v.matches_reference = matches_reference
    v
end

@testset "classify_allele distinguishes reference/snp/indel" begin
    @test classify_allele("A", "A")     == :reference
    @test classify_allele("A", "G")     == :snp
    @test classify_allele("ACA", "A")   == :indel   # deletion
    @test classify_allele("A", "AT")    == :indel   # insertion
end

@testset "aggregate_locus_alleles keys by (ref,base), no collapse" begin
    vars = [
        mkvar(strain="S1", reference="A",   base="G"),
        mkvar(strain="S2", reference="ACA", base="A"),
        mkvar(strain="REF", reference="A",  base="A", matches_reference=1),
    ]
    (stats, total) = aggregate_locus_alleles(vars)
    @test total == 3
    @test haskey(stats, ("ACA", "A"))
    @test haskey(stats, ("A", "A"))
    @test haskey(stats, ("A", "G"))
    @test stats[("ACA","A")].weight == 1
    @test length(stats[("A","G")].strains) == 1
end

@testset "aggregate: complex decomposition counts each strain's ploidy once" begin
    # LmjF.01:85879 shape — 2 diploid strains carry BOTH a SNP (T>C) and a
    # deletion (TGT>T) via two 1/1 records; 2 ref strains + 1 haploid synthetic ref.
    mkslots(strain, ref, base, ploidy) = begin
        v = Variation(); v.strain = strain; v.reference = ref; v.base = base
        v.ploidy = ploidy; v.percent = "100"; v.coverage = "10"
        v.allele_slots = fill(base, ploidy); v
    end
    vars = [
        mkslots("LV39",     "T",   "C", 2), mkslots("LV39",     "TGT", "T", 2),
        mkslots("LV39cl5",  "T",   "C", 2), mkslots("LV39cl5",  "TGT", "T", 2),
        mkslots("Fried",    "T",   "T", 2), mkslots("Seid",     "T",   "T", 2),
        mkslots("synthref", "T",   "T", 1),
    ]
    (stats, total) = aggregate_locus_alleles(vars)
    @test total == 9                         # 2+2+2+2+1, each strain once
    @test stats[("T","C")].weight   == 4     # both chromosomes of 2 strains
    @test stats[("TGT","T")].weight == 4
    @test stats[("T","T")].weight   == 5     # Fried 2 + Seid 2 + synthref 1
end

@testset "aggregate: split 1/2 compound het fabricates no reference" begin
    # LmjF.01:8962 shape — Seidman is 1/2 (TA/TAA), split into 1/. and ./1.
    het(strain, ref, alt) = begin
        v = Variation(); v.strain = strain; v.reference = ref; v.base = alt
        v.ploidy = 2; v.percent = "100"; v.coverage = "12"
        v.allele_slots = [alt]; v          # one non-missing slot, no ref
    end
    refstrain(strain, ploidy) = begin
        v = Variation(); v.strain = strain; v.reference = "T"; v.base = "T"
        v.ploidy = ploidy; v.percent = "100"; v.coverage = "20"
        v.allele_slots = fill("T", ploidy); v
    end
    vars = [
        het("Seid", "T", "TA"), het("Seid", "T", "TAA"),
        refstrain("LV39", 2), refstrain("Fried", 2), refstrain("LV39cl5", 2),
        refstrain("synthref", 1),
    ]
    (stats, total) = aggregate_locus_alleles(vars)
    @test total == 9                          # Seid counted once as ploidy 2
    @test stats[("T","TA")].weight  == 1
    @test stats[("T","TAA")].weight == 1
    @test stats[("T","T")].weight   == 7      # Seid contributes 0 reference
    @test !("Seid" in stats[("T","T")].strains)
end

@testset "chromosome_alleles: legacy hom derives from base×ploidy" begin
    v = Variation(); v.reference = "A"; v.base = "G"; v.ploidy = 2
    @test chromosome_alleles(v) == ["G", "G"]
end

@testset "chromosome_alleles: legacy het derives [ref, alt]" begin
    v = Variation(); v.reference = "A"; v.base = "R"; v.alt_allele = "G"; v.ploidy = 2
    @test chromosome_alleles(v) == ["A", "G"]
end

@testset "chromosome_alleles: explicit allele_slots wins over legacy fields" begin
    v = Variation(); v.reference = "T"; v.base = "TA"; v.ploidy = 2
    v.allele_slots = ["TA"]
    @test chromosome_alleles(v) == ["TA"]
end

@testset "gt_to_base: half-missing 1/. returns the present alt" begin
    @test gt_to_base("1/.", "T", ["TA"]) == "TA"
    @test gt_to_base("./1", "T", ["TA"]) == "TA"
end

@testset "nonref_alt_alleles: half-missing keeps the present alt" begin
    @test nonref_alt_alleles("1/.", ["TA"]) == ["TA"]
    @test nonref_alt_alleles("./1", ["TA"]) == ["TA"]
end

@testset "compute_percent: half-missing uses present alt AO" begin
    fmt = Dict("AO" => "7", "RO" => "0")
    @test compute_percent(fmt, "1/.") == "100.00"
    @test compute_percent(fmt, "./1") == "100.00"
end

@testset "build_variations_from_record: 1/. yields one alt slot, no ref" begin
    rec = make_vcf_record(pos=8962, ref="T", alts=["TA"],
                          format_keys=["GT","AO","RO"], sample_data=["1/.:7:0"])
    cov = make_coverage("s1", 8961, 200, 12.0)
    vars = build_variations_from_record(rec, ["s1"], cov, 2)
    @test length(vars) == 1
    @test vars[1].allele_slots == ["TA"]
    @test vars[1].ploidy == 2
end

# ---------------------------------------------------------------------------
# n-ploidy genotypes (e.g. triploid GTs from aneuploid tritryp/fungal strains,
# preserved by bcftools merge). See merged VCF Chr1_A_fumigatus_Af293:12263,
# which mixes diploid 1/1 with triploid ././., 1/./., ./1/., ././1 in one record.
# ---------------------------------------------------------------------------

@testset "is_missing_gt: true only when every slot is missing, any ploidy" begin
    @test is_missing_gt("")       == true
    @test is_missing_gt(".")      == true
    @test is_missing_gt("./.")    == true
    @test is_missing_gt(".|.")    == true
    @test is_missing_gt("././.")  == true
    @test is_missing_gt(".|.|.")  == true
    @test is_missing_gt("1/./.")  == false
    @test is_missing_gt("./1/.")  == false
    @test is_missing_gt("0/0/0")  == false
    @test is_missing_gt("1/1")    == false
end

@testset "gt_to_base: triploid all-missing is empty" begin
    @test gt_to_base("././.", "T", ["TA"]) == ""
    @test gt_to_base(".|.|.", "T", ["TA"]) == ""
end

@testset "gt_to_base: triploid with one present allele returns that allele" begin
    @test gt_to_base("1/./.", "T", ["TA"]) == "TA"
    @test gt_to_base("./1/.", "T", ["TA"]) == "TA"
    @test gt_to_base("././1", "T", ["TA"]) == "TA"
end

@testset "gt_to_base: homozygous triploid returns the single allele" begin
    @test gt_to_base("1/1/1", "T", ["TA"]) == "TA"   # hom-alt
    @test gt_to_base("0/0/0", "A", ["G"])  == "A"    # hom-ref
end

@testset "gt_to_base: triploid het of two single-char SNP alleles → IUPAC" begin
    @test gt_to_base("0/0/1", "A", ["G"]) == "R"     # {A,G} → R, ref-heavy
    @test gt_to_base("0/1/1", "A", ["G"]) == "R"     # {A,G} → R, alt-heavy
end

@testset "nonref_alt_alleles: n-ploidy" begin
    @test nonref_alt_alleles("././.", ["TA"]) == String[]
    @test nonref_alt_alleles("1/./.", ["TA"]) == ["TA"]
    @test nonref_alt_alleles("1/1/1", ["TA"]) == ["TA"]
    @test nonref_alt_alleles("1/2/.", ["TA","TAA"]) == ["TA","TAA"]
end

@testset "compute_percent: n-ploidy" begin
    fmt = Dict("AO" => "68", "RO" => "0")
    @test compute_percent(fmt, "1/./.") == "100.00"
    @test compute_percent(fmt, "././.") == "0.0"
    @test compute_percent(fmt, "0/0/0") == "0.0"
end

@testset "build_variations_from_record: triploid 1/./. yields one alt slot, ploidy 3" begin
    # Real-data shape: GT:DP:AD:RO:QR:AO:QA with a triploid partial-missing call
    rec = make_vcf_record(pos=12263, ref="G", alts=["A"],
                          format_keys=["GT","DP","AD","RO","QR","AO","QA"],
                          sample_data=["1/./.:128:0,68:0:0:68:2529"])
    cov = make_coverage("s1", 12262, 200, 128.0)
    vars = build_variations_from_record(rec, ["s1"], cov, 2)
    @test length(vars) == 1
    @test vars[1].allele_slots == ["A"]
    @test vars[1].ploidy == 3
end

@testset "build_variations_from_record: triploid all-missing does not crash, synthesizes ref when covered" begin
    rec = make_vcf_record(pos=12263, ref="G", alts=["A"],
                          format_keys=["GT","DP","AD","RO","QR","AO","QA"],
                          sample_data=["././.:.:.:.:.:.:."])
    cov = make_coverage("s1", 12000, 13000, 40.0)
    vars = build_variations_from_record(rec, ["s1"], cov, 3)
    @test length(vars) == 1
    @test vars[1].matches_reference == 1
    @test vars[1].base == "G"
    @test vars[1].ploidy == 3           # from --ploidy default, no GT to read
end

@testset "build_variations_from_record: homozygous triploid is not treated as het" begin
    # 1/1/1 is hom-alt, not het — must behave like diploid 1/1 (alt_allele empty).
    rec = make_vcf_record(pos=12263, ref="G", alts=["A"],
                          format_keys=["GT","DP","AD","RO","QR","AO","QA"],
                          sample_data=["1/1/1:90:0,90:0:0:90:3000"])
    cov = make_coverage("s1", 12000, 13000, 90.0)
    vars = build_variations_from_record(rec, ["s1"], cov, 2)
    @test length(vars) == 1
    @test vars[1].base == "A"
    @test vars[1].ploidy == 3
    @test vars[1].allele_slots == ["A","A","A"]
    @test vars[1].alt_allele == ""          # hom-alt ⇒ no het alt component
end

@testset "gt_to_ca: emits one CA slot per GT slot at any ploidy" begin
    # diploid baseline (locks existing behavior)
    @test gt_to_ca("0/1", 1, ["k0"], ["r0"]) == "r0/k0"
    @test gt_to_ca("./.", 1, ["k0"], ["r0"]) == "."          # missing → single dot
    @test gt_to_ca("0|1", 1, ["k0"], ["r0"]) == "r0|k0"      # phased separator preserved
    # triploid: every slot must appear, separators preserved, '.' slot stays '.'
    @test gt_to_ca("0/1/.", 1, ["k0"], ["r0"]) == "r0/k0/."
    @test gt_to_ca("1/1/1", 1, ["k0"], ["r0"]) == "k0/k0/k0"
    @test gt_to_ca("0/0/0", 1, ["k0"], ["r0"]) == "r0/r0/r0"
end

@testset "write_snp_feature emits per-class SNP+indel columns without collapse" begin
    vars = [
        mkvar(strain="S1", reference="A",   base="G", coverage="30", percent="100"),
        mkvar(strain="S2", reference="ACA", base="A", coverage="20", percent="100"),
        mkvar(strain="REF", reference="A",  base="A", coverage="0",  percent="100", matches_reference=1),
    ]
    buf = IOBuffer()
    write_snp_feature(buf, vars, 1, "LmjF.01", 13850, "REF", ["S1","S2"])
    cols = split(chomp(String(take!(buf))), '\t')
    @test length(cols) == 31
    @test cols[1] == "13850"
    @test cols[5] == "MIXED"
    @test cols[13] == "A"                       # snp_ref_allele
    @test cols[14] == "A"                       # snp_major_allele = reference (tie, "A"<"G")
    @test cols[20] == ""                        # snp_major_genomic_hgvs blank for ref slot
    @test cols[17] == "G"                       # snp_minor_allele = the SNP alt
    @test cols[21] == "LmjF.01:g.13850A>G"      # snp_minor_genomic_hgvs
    @test cols[22] == "ACA"                     # indel_ref_allele
    @test cols[23] == "A"                       # indel_major_allele = deletion (no indel ref key present)
    @test occursin("del", cols[29])             # indel_major_genomic_hgvs
    @test cols[31] == "frameshift"              # indel_frame_effect (1-3=-2)
end

@testset "write_snp_feature SNP-only locus leaves indel family empty" begin
    vars = [
        mkvar(strain="S1", reference="A", base="G"),
        mkvar(strain="REF", reference="A", base="A", matches_reference=1),
    ]
    buf = IOBuffer()
    write_snp_feature(buf, vars, 1, "chr1", 100, "REF", ["S1"])
    cols = split(chomp(String(take!(buf))), '\t')
    @test cols[5] == "SNV"
    @test cols[22] == ""                        # indel_ref_allele empty
    @test cols[31] == ""                        # indel_frame_effect empty
end

@testset "write_allele_file keeps deletion and reference as distinct rows" begin
    vars = [
        mkvar(strain="S1", reference="A",   base="G"),
        mkvar(strain="S2", reference="ACA", base="A"),
        mkvar(strain="REF", reference="A",  base="A", matches_reference=1),
    ]
    buf = IOBuffer()
    sid = Dict("S1"=>1, "S2"=>2, "REF"=>3)
    write_allele_file(buf, vars, "LmjF.01", 13850, sid)
    rows = [split(l, '\t') for l in filter(!isempty, split(String(take!(buf)), "\n"))]
    del  = [r for r in rows if r[3] == "A" && r[9] == "0"]
    refr = [r for r in rows if r[3] == "A" && r[9] == "1"]
    @test length(del)  == 1                 # deletion A, matches_reference 0
    @test length(refr) == 1                 # reference A, matches_reference 1
    @test occursin("del", del[1][10])       # deletion has a del g.HGVS
    @test refr[1][10] == "."                # reference row g.HGVS is "."
end

@testset "remap_sample_for_split marks other-alt slots missing, not reference" begin
    # n_orig_alts=2, target_alt_i=2: slot "1" is a NON-target alt → "." (its call
    # lives in the sibling split record); "." stays "."
    @test remap_sample_for_split("1/.", ["GT"], 2, 2) == "./."
    @test remap_sample_for_split("./1", ["GT"], 2, 2) == "./."
    # target_alt_i=1: slot "1" IS the target alt → 1; "." stays "."
    @test remap_sample_for_split("1/.", ["GT"], 2, 1) == "1/."
    # reference slot stays reference
    @test remap_sample_for_split("0/1", ["GT"], 2, 1) == "0/1"
    @test remap_sample_for_split("0/1", ["GT"], 2, 2) == "0/."
    # full-missing GT is left as-is (guarded earlier)
    @test remap_sample_for_split("./.", ["GT"], 2, 2) == "./."
    # sanity: existing biallelic behavior unaffected (n_orig_alts=1 → unchanged)
    @test remap_sample_for_split("1/2", ["GT"], 1, 1) == "1/2"
end

@testset "remap_sample_for_split remaps triploid GTs, other-alt → missing" begin
    # Triploid split multiallelic. Target alt → 1, reference (0) stays 0,
    # every other alt index → "." ; '.' slots and both separators preserved.
    @test remap_sample_for_split("./2/.", ["GT"], 3, 2) == "./1/."   # target alt → 1
    @test remap_sample_for_split("./2/.", ["GT"], 3, 1) == "././."   # non-target alt → .
    @test remap_sample_for_split("1/./.", ["GT"], 3, 1) == "1/./."   # target alt → 1
    @test remap_sample_for_split("1/2/3", ["GT"], 3, 2) == "./1/."   # only target survives
    @test remap_sample_for_split("0/2/0", ["GT"], 3, 2) == "0/1/0"   # ref slots stay ref
    @test remap_sample_for_split("././.", ["GT"], 3, 2) == "././."   # full-missing untouched
end

@testset "write_snp_feature het_strain_count counts distinct strains" begin
    # One strain, two het records (split 1/2). het_strain_count must be 1, not 2.
    v1 = Variation(); v1.strain="Seid"; v1.reference="T"; v1.base="TA";  v1.alt_allele="TA";  v1.ploidy=2; v1.coverage="12"; v1.percent="58"
    v2 = Variation(); v2.strain="Seid"; v2.reference="T"; v2.base="TAA"; v2.alt_allele="TAA"; v2.ploidy=2; v2.coverage="12"; v2.percent="42"
    buf = IOBuffer()
    write_snp_feature(buf, [v1, v2], 0, "LmjF.01", 8962, "synthref", ["Seid"])
    fields = split(strip(String(take!(buf))), "\t")
    @test fields[12] == "1"   # het_strain_count (column 12)
end
