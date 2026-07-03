#!/usr/bin/env julia

# processSequenceVariations.jl
# Reads a merged multi-sample FreeBayes VCF and a coordinate-sorted TSV cache file,
# streams them concurrently in a sorted merge, annotates coding variants via SQLite
# transcript/indel databases, and writes five output files:
#   output.vcf (CANN-annotated VCF for SnpEff), cache.tsv (position cache),
#   snpFeature.dat, allele.dat, product.dat
# Coverage information is read from a coverage.tsv produced by mergeCoverageBeds.

using SQLite
using SQLite.DBInterface: execute
using Printf

# ---------------------------------------------------------------------------
# Global debug flag
# ---------------------------------------------------------------------------

global DEBUG = false

function debug_log(msg...)
    DEBUG && println(stderr, "[DEBUG] ", msg...)
end

# ---------------------------------------------------------------------------
# CLI argument parsing (hand-rolled, no package dependency)
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Standard codon table (NCBI translation table 1)
# ---------------------------------------------------------------------------

const CODON_TABLE = Dict{String,String}(
    "TTT"=>"F","TTC"=>"F","TTA"=>"L","TTG"=>"L",
    "CTT"=>"L","CTC"=>"L","CTA"=>"L","CTG"=>"L",
    "ATT"=>"I","ATC"=>"I","ATA"=>"I","ATG"=>"M",
    "GTT"=>"V","GTC"=>"V","GTA"=>"V","GTG"=>"V",
    "TCT"=>"S","TCC"=>"S","TCA"=>"S","TCG"=>"S",
    "CCT"=>"P","CCC"=>"P","CCA"=>"P","CCG"=>"P",
    "ACT"=>"T","ACC"=>"T","ACA"=>"T","ACG"=>"T",
    "GCT"=>"A","GCC"=>"A","GCA"=>"A","GCG"=>"A",
    "TAT"=>"Y","TAC"=>"Y","TAA"=>"*","TAG"=>"*",
    "CAT"=>"H","CAC"=>"H","CAA"=>"Q","CAG"=>"Q",
    "AAT"=>"N","AAC"=>"N","AAA"=>"K","AAG"=>"K",
    "GAT"=>"D","GAC"=>"D","GAA"=>"E","GAG"=>"E",
    "TGT"=>"C","TGC"=>"C","TGA"=>"*","TGG"=>"W",
    "CGT"=>"R","CGC"=>"R","CGA"=>"R","CGG"=>"R",
    "AGT"=>"S","AGC"=>"S","AGA"=>"R","AGG"=>"R",
    "GGT"=>"G","GGC"=>"G","GGA"=>"G","GGG"=>"G"
)

# ---------------------------------------------------------------------------
# IUPAC ambiguity tables
# ---------------------------------------------------------------------------

const IUPAC_EXPAND = Dict{Char,Vector{Char}}(
    'A'=>['A'], 'G'=>['G'], 'C'=>['C'], 'T'=>['T'],
    'R'=>['A','G'], 'Y'=>['C','T'], 'K'=>['G','T'], 'M'=>['A','C'],
    'S'=>['G','C'], 'W'=>['A','T'],
    'B'=>['G','T','C'], 'D'=>['G','A','T'], 'H'=>['A','C','T'], 'V'=>['G','C','A'],
    'N'=>['A','G','C','T']
)

# Reverse IUPAC: Set of two bases -> ambiguity code (for diploid het calls)
const IUPAC_COMPRESS = Dict{Set{Char},Char}(
    Set(['A','G'])=>'R', Set(['C','T'])=>'Y', Set(['G','T'])=>'K',
    Set(['A','C'])=>'M', Set(['G','C'])=>'S', Set(['A','T'])=>'W',
)

# ---------------------------------------------------------------------------
# HSSS (HighSpeedSnpSearch) binary strain file support
# ---------------------------------------------------------------------------

const HSSS_CUTOFFS = (20, 40, 60, 80)

const HSSS_ALLELE_CODE = Dict{Char,Int8}(
    'A'=>Int8(1), 'a'=>Int8(1),
    'C'=>Int8(2), 'c'=>Int8(2),
    'G'=>Int8(3), 'g'=>Int8(3),
    'T'=>Int8(4), 't'=>Int8(4),
)

mutable struct HsssState
    ref_fhs::Vector{IO}               # one referenceGenome.dat handle per cutoff
    contig_fhs::Vector{IO}            # one contigIdToSourceId.dat handle per cutoff
    strain_fhs::Vector{Dict{Int,IO}}  # per cutoff: strain_index -> file handle (non-ref only)
    strain_index::Dict{String,Int}    # strain_name -> integer index (ref=1, others 2..N)
    seq_index::Int
    current_seq_id::String
end

function open_hsss_writers(reference_strain::String, all_strains::Vector{String},
                            base_dir::String=".")::HsssState
    non_ref = filter(s -> s != reference_strain, all_strains)
    strain_index = Dict{String,Int}(reference_strain => 1)
    for (i, s) in enumerate(non_ref)
        strain_index[s] = i + 1
    end

    ref_fhs    = IO[]
    contig_fhs = IO[]
    strain_fhs = Dict{Int,IO}[]

    for cutoff in HSSS_CUTOFFS
        dir = joinpath(base_dir, "hsss_readFreq$(cutoff)")
        mkpath(dir)

        open(joinpath(dir, "strainIdToName.dat"), "w") do f
            write(f, "1\t$(reference_strain)\n")
            for (i, s) in enumerate(non_ref)
                write(f, "$(i+1)\t$(s)\n")
            end
        end

        push!(ref_fhs,    open(joinpath(dir, "referenceGenome.dat"), "w"))
        push!(contig_fhs, open(joinpath(dir, "contigIdToSourceId.dat"), "w"))

        sfhs = Dict{Int,IO}()
        for (i, s) in enumerate(non_ref)
            sfhs[i + 1] = open(joinpath(dir, string(i + 1)), "w")
        end
        push!(strain_fhs, sfhs)
    end

    HsssState(ref_fhs, contig_fhs, strain_fhs, strain_index, 0, "")
end

function close_hsss_writers(state::HsssState)
    for fh in state.ref_fhs;    close(fh); end
    for fh in state.contig_fhs; close(fh); end
    for sfhs in state.strain_fhs
        for (_, fh) in sfhs; close(fh); end
    end
end

function write_hsss_record(fh::IO, seq_idx::Int, location::Int, allele_c::Int8, product_c::Int8)
    write(fh, htol(Int16(seq_idx)))
    write(fh, htol(Int32(location)))
    write(fh, allele_c)
    write(fh, product_c)
end

# write_hsss_position! is defined after the Variation struct (see below)

"""
    expand_codon(codon::String) -> Vector{String}

Cartesian product expansion of IUPAC ambiguity codes in a 3-character codon.
"""
function expand_codon(codon::String)
    c = uppercase(codon)
    if length(c) != 3
        c = "NNN"
    end
    bases = [get(IUPAC_EXPAND, ch, ['N']) for ch in c]
    codons = String[]
    for b1 in bases[1], b2 in bases[2], b3 in bases[3]
        push!(codons, string(b1, b2, b3))
    end
    codons
end

function translate_codon(codon::String)
    get(CODON_TABLE, uppercase(codon), "X")
end

# ---------------------------------------------------------------------------
# GTF parsing → CDS interval structures
# ---------------------------------------------------------------------------

struct CDSInterval
    seq_id::String
    start_pos::Int
    end_pos::Int
    strand::Int
    transcript_id::String
    cds_number::Int
end

struct TranscriptInfo
    seq_id::String
    strand::Int
    exons::Vector{CDSInterval}   # sorted by start_pos
end

function parse_gtf(gtf_file::String)
    debug_log("Parsing GTF file: ", gtf_file)
    transcript_cds = Dict{String, Vector{Tuple{String,Int,Int,Int}}}()
    transcript_order = String[]

    open(gtf_file, "r") do fh
        for line in eachline(fh)
            startswith(line, '#') && continue
            fields = split(line, '\t')
            length(fields) < 9 && continue
            fields[3] == "CDS" || continue

            seq_id = fields[1]
            start_pos = parse(Int, fields[4])
            end_pos = parse(Int, fields[5])
            strand = (fields[7] == "-") ? -1 : 1
            attributes = fields[9]

            tid = ""
            for attr in split(attributes, ';')
                attr = strip(attr)
                if startswith(attr, "transcript_id")
                    parts = split(attr, ' '; limit=2)
                    if length(parts) == 2
                        tid = replace(parts[2], "\"" => "")
                    end
                    break
                end
            end
            isempty(tid) && continue

            if !haskey(transcript_cds, tid)
                transcript_cds[tid] = Tuple{String,Int,Int,Int}[]
                push!(transcript_order, tid)
            end
            push!(transcript_cds[tid], (seq_id, start_pos, end_pos, strand))
        end
    end

    all_intervals = CDSInterval[]
    transcript_info = Dict{String, TranscriptInfo}()

    for tid in transcript_order
        entries = transcript_cds[tid]
        sort!(entries; by = e -> e[2])

        seq_id = entries[1][1]
        strand = entries[1][4]

        exons = CDSInterval[]
        for (idx, (sid, s, e, st)) in enumerate(entries)
            interval = CDSInterval(sid, s, e, st, tid, idx)
            push!(exons, interval)
            push!(all_intervals, interval)
        end

        transcript_info[tid] = TranscriptInfo(seq_id, strand, exons)
    end

    sort!(all_intervals; by = iv -> (iv.seq_id, iv.start_pos))
    debug_log("GTF parsed: ", length(transcript_info), " transcripts, ",
              length(all_intervals), " CDS intervals")
    (all_intervals, transcript_info)
end

"""
    find_all_cds(intervals, seq_id, location) -> Vector{CDSInterval}

Returns all CDS intervals on seq_id that overlap location.
intervals must be sorted by (seq_id, start_pos) ascending.
"""
function find_all_cds(intervals::Vector{CDSInterval}, seq_id::String, location::Int)
    # Binary search: rightmost interval with start_pos <= location on this seq_id
    lo, hi = 1, length(intervals)
    candidate = 0
    while lo <= hi
        mid = (lo + hi) ÷ 2
        iv = intervals[mid]
        if iv.seq_id < seq_id || (iv.seq_id == seq_id && iv.start_pos <= location)
            candidate = mid
            lo = mid + 1
        else
            hi = mid - 1
        end
    end

    candidate == 0 && return CDSInterval[]

    # Scan backward from candidate collecting all intervals where end_pos >= location.
    # All have start_pos <= location by construction; we just need end_pos >= location.
    hits = CDSInterval[]
    i = candidate
    while i >= 1
        iv = intervals[i]
        iv.seq_id != seq_id && break
        iv.end_pos >= location && push!(hits, iv)
        i -= 1
    end
    hits
end

function compute_position_in_cds(tinfo::TranscriptInfo, location::Int, cds_number::Int)
    prior_len = 0
    if tinfo.strand == 1
        for exon in tinfo.exons
            if exon.cds_number == cds_number
                return prior_len + (location - exon.start_pos) + 1
            end
            prior_len += exon.end_pos - exon.start_pos + 1
        end
    else
        # Negative strand: CDS order is highest→lowest genomic coordinate.
        # Iterate exons in reverse genomic order; within-exon offset measures
        # distance from the 3' end of the exon on the genome (= 5' direction in CDS).
        for i in length(tinfo.exons):-1:1
            exon = tinfo.exons[i]
            if exon.cds_number == cds_number
                return prior_len + (exon.end_pos - location) + 1
            end
            prior_len += exon.end_pos - exon.start_pos + 1
        end
    end
    error("cds_number $cds_number out of range for transcript")
end

function position_in_codon(pos_in_cds::Int)
    r = pos_in_cds % 3
    r == 0 ? 3 : r
end

function extract_codon(sequence::String, pos_in_cds::Int)
    pic = position_in_codon(pos_in_cds)
    codon_start = pos_in_cds - pic  # 0-based index of first base
    if codon_start + 3 > length(sequence)
        return "NNN"
    end
    sequence[codon_start+1 : codon_start+3]
end

# ---------------------------------------------------------------------------
# Frameshift precomputation
# ---------------------------------------------------------------------------

function precompute_frameshifts(indel_db::SQLite.DB, transcript_info::Dict{String,TranscriptInfo})
    debug_log("Precomputing frameshift information...")
    fs_info = Dict{String, Dict{String, Tuple{Bool, Int}}}()

    current_strain = ""
    current_tid = ""
    cumsum = 0
    found_frameshift = false
    frameshift_pos = 0

    function flush_current()
        if !isempty(current_strain) && !isempty(current_tid)
            if !haskey(fs_info, current_strain)
                fs_info[current_strain] = Dict{String, Tuple{Bool, Int}}()
            end
            fs_info[current_strain][current_tid] = (found_frameshift, frameshift_pos)
        end
    end

    for row in execute(indel_db,
        "SELECT strain, transcript_id, position, shift_amount FROM indels ORDER BY strain, transcript_id, position")
        strain = row[1]::String
        tid = row[2]::String
        pos = row[3]::Int
        shift = row[4]::Int

        if strain != current_strain || tid != current_tid
            flush_current()
            current_strain = strain
            current_tid = tid
            cumsum = 0
            found_frameshift = false
            frameshift_pos = 0
        end

        cumsum += shift
        if !found_frameshift && (cumsum % 3 != 0)
            found_frameshift = true
            frameshift_pos = pos
        end
    end
    flush_current()

    fs_count = sum(length(v) for v in values(fs_info); init=0)
    debug_log("Frameshift precomputation complete: ", fs_count, " strain-transcript pairs with indels")
    fs_info
end

function check_downstream_of_frameshift(
    fs_info::Dict{String, Dict{String, Tuple{Bool, Int}}},
    strain::String,
    transcript::String,
    location::Int
)
    strain_fs = get(fs_info, strain, nothing)
    isnothing(strain_fs) && return 0
    tid_fs = get(strain_fs, transcript, nothing)
    isnothing(tid_fs) && return 0
    has_fs, fs_pos = tid_fs
    has_fs || return 0
    location > fs_pos ? 1 : 0
end

# ---------------------------------------------------------------------------
# SQLite sequence loading
# ---------------------------------------------------------------------------

function load_transcript_sequences(db::SQLite.DB, transcript_id::String)
    result = Dict{String,String}()
    for row in execute(db, "SELECT strain, sequence FROM coding_sequences WHERE transcript_id = ?", [transcript_id])
        result[row[1]::String] = row[2]::String
    end
    result
end

function get_indel_shift(db::SQLite.DB, transcript_id::String, strain::String, position::Int)
    row = first(execute(db,
        "SELECT COALESCE(SUM(shift_amount), 0) FROM indels WHERE transcript_id = ? AND strain = ? AND position < ?",
        [transcript_id, strain, position]))
    row[1]::Int
end

# ---------------------------------------------------------------------------
# VCF record structure
# ---------------------------------------------------------------------------

struct VCFRecord
    chrom::String
    pos::Int
    ref::String
    alts::Vector{String}
    info::String
    format_keys::Vector{String}
    sample_data::Vector{String}   # raw per-sample FORMAT strings
end

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Variation record (mutable, used during processing)
# ---------------------------------------------------------------------------

mutable struct Variation
    sequence_source_id::String
    location::Int
    strain::String
    reference::String
    base::String
    alt_allele::String   # actual VCF alt base for het calls; empty otherwise
    coverage::String
    percent::String      # alt_fraction = AO/(RO+AO)*100; for het, ref_fraction = 100 - percent
    quality::String
    pvalue::String
    snp_source_id::String
    is_coding::Int
    position_in_cds::Int
    position_in_codon::Int
    downstream_of_frameshift::Int
    transcript::String
    product::Vector{String}
    reference_codon::String
    codon::String
    has_nonsynonomous::Int
    cds_number::Int
    matches_reference::Int
    ploidy::Int
end

function Variation()
    Variation("", 0, "", "", "", "", "", "", "", "", "",
              0, 0, 0, 0, "", String[], "", "", 0, 0, 0, 1)
end

# ---------------------------------------------------------------------------
# HSSS write_hsss_position! — defined here because it references Variation
# ---------------------------------------------------------------------------

# Same length-comparison criterion as has_snp_allele, applied to an already-built Variation
# rather than a raw VCF record.
is_snp_variation(v::Variation) = length(v.reference) == length(v.base)

function write_hsss_position!(
    state::HsssState,
    variations::Vector{Variation},
    reference_strain::String,
    seq_id::String,
    location::Int,
    all_strains::Vector{String},
    product_code::Int8
)
    # Update sequence index when chromosome changes
    if seq_id != state.current_seq_id
        state.seq_index     += 1
        state.current_seq_id = seq_id
        for fh in state.contig_fhs
            write(fh, "$(state.seq_index)\t$(seq_id)\n")
        end
    end
    seq_idx = state.seq_index

    # Index variations by strain
    found = Dict{String, Vector{Variation}}()
    for v in variations
        push!(get!(found, v.strain, Variation[]), v)
    end

    ref_vars     = get(found, reference_strain, Variation[])
    ref_base     = isempty(ref_vars) ? "" : ref_vars[1].base
    ref_allele_c = get(HSSS_ALLELE_CODE, isempty(ref_base) ? ' ' : ref_base[1], Int8(0))

    for (ci, cutoff) in enumerate(HSSS_CUTOFFS)
        # Determine which non-ref strains pass the cutoff
        passing = Set{String}()
        for (strain, svars) in found
            strain == reference_strain && continue
            if any(v -> !isempty(v.percent) && parse(Float64, v.percent) >= cutoff, svars)
                push!(passing, strain)
            end
        end

        # Skip position for this cutoff if no passing strain differs from reference
        has_notable = any(passing) do strain
            any(v -> v.matches_reference != 1, get(found, strain, Variation[]))
        end
        has_notable || continue

        # Write reference genome record
        write_hsss_record(state.ref_fhs[ci], seq_idx, location, ref_allele_c, product_code)

        # Write per-strain records
        for strain in all_strains
            strain == reference_strain && continue
            sidx = get(state.strain_index, strain, 0)
            sidx == 0 && continue
            sfh = get(state.strain_fhs[ci], sidx, nothing)
            isnothing(sfh) && continue

            svars = get(found, strain, nothing)

            if isnothing(svars)
                # Strain absent from this position
                write_hsss_record(sfh, seq_idx, location, Int8(0), Int8(0))
                continue
            end

            if strain ∉ passing
                # Present but below cutoff: treat as unknown
                write_hsss_record(sfh, seq_idx, location, Int8(0), Int8(0))
                continue
            end

            # HSSS is SNP-only: drop indel-length alleles rather than writing a bogus
            # allele code of 0 for them (see has_snp_allele/is_snp_variation).
            snp_svars = filter(is_snp_variation, svars)
            if isempty(snp_svars)
                # Strain's only allele(s) at this position are indels: treat as unknown
                write_hsss_record(sfh, seq_idx, location, Int8(0), Int8(0))
                continue
            end

            # Het strains write one record per allele (including ref-matching); matches Perl hsssCreateStrainFiles behavior
            is_het = length(snp_svars) > 1
            for sv in snp_svars
                # Skip non-het, reference-matching variations
                !is_het && sv.matches_reference == 1 && continue
                allele_c = get(HSSS_ALLELE_CODE, isempty(sv.base) ? ' ' : sv.base[1], Int8(0))
                write_hsss_record(sfh, seq_idx, location, allele_c, product_code)
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Processing data structures
# ---------------------------------------------------------------------------

struct ProcessingContext
    reference_strain::String
    undone_strains::Set{String}
    cds_intervals::Vector{CDSInterval}
    transcript_info::Dict{String,TranscriptInfo}
    transcript_db::SQLite.DB
    indel_db::SQLite.DB
    fs_info::Dict{String, Dict{String, Tuple{Bool, Int}}}
    all_strains::Vector{String}
    sample_id_map::Dict{String,Int}
    ploidy::Int
end

struct OutputWriters
    vcf_fh::IO
    snp_fh::IO
    allele_fh::IO
    transcript_product_fh::IO
    hsss::HsssState
end

function write_sample_dat(all_strains::Vector{String}, path::String="sample.dat")
    open(path, "w") do fh
        write(fh, "strain_id\tsample_name\n")
        for (i, name) in enumerate(all_strains)
            write(fh, "$(i)\t$(name)\n")
        end
    end
end

struct PositionAnnotation
    is_coding::Int
    transcript_id::String
    cds_number::Int
    pos_in_cds::Int
    pos_in_codon_val::Int
    ref_codon::String
    ref_product::String
end

mutable struct TranscriptSequenceCache
    seqs::Dict{String, Dict{String,String}}   # transcript_id -> strain -> sequence
end

function load_transcript_if_needed!(cache::TranscriptSequenceCache, db::SQLite.DB, transcript_id::String)
    haskey(cache.seqs, transcript_id) && return
    cache.seqs[transcript_id] = load_transcript_sequences(db, transcript_id)
end

function get_strain_seq(cache::TranscriptSequenceCache, transcript_id::String, strain::String)::String
    t = get(cache.seqs, transcript_id, nothing)
    isnothing(t) ? "" : get(t, strain, "")
end

# ---------------------------------------------------------------------------
# PeekedFile: one-line lookahead over any IO (file, subprocess pipe, IOBuffer)
# ---------------------------------------------------------------------------

mutable struct PeekedFile
    fh::IO
    line::String
    exhausted::Bool
end

function advance!(pf::PeekedFile)
    if eof(pf.fh)
        pf.line = ""
        pf.exhausted = true
    else
        pf.line = readline(pf.fh)
        if isempty(pf.line) && eof(pf.fh)
            pf.exhausted = true
        end
    end
end

function close_peeked(pf::PeekedFile)
    close(pf.fh)
end

# ---------------------------------------------------------------------------
# VCF I/O
# ---------------------------------------------------------------------------

"""
    parse_vcf_header(io) -> (all_strains, chrom_rank, info_headers)

Reads ## meta lines, builds chrom_rank from ##contig lines, extracts sample
names from #CHROM line. Leaves io positioned at first data line.
"""
function parse_vcf_header(io::IO)
    chrom_rank   = Dict{String,Int}()
    all_strains  = String[]
    info_headers = String[]
    contig_count = 0

    for line in eachline(io)
        if startswith(line, "##INFO")
            push!(info_headers, line)
        elseif startswith(line, "##contig")
            m = match(r"##contig=<ID=([^,>]+)", line)
            if !isnothing(m)
                contig_count += 1
                chrom_rank[m.captures[1]] = contig_count
            end
        elseif startswith(line, "#CHROM")
            fields = split(line, '\t')
            all_strains = String[String(fields[i]) for i in 10:length(fields)]
            break
        end
    end

    debug_log("VCF header: ", length(all_strains), " samples, ",
              length(chrom_rank), " contigs")
    (all_strains, chrom_rank, info_headers)
end

"""
    open_vcf_peeked(path) -> (PeekedFile, all_strains, chrom_rank, info_headers)

Opens a bgzip-compressed VCF via subprocess, parses its header, returns
a PeekedFile positioned at the first data line.
"""
function open_vcf_peeked(path::String)
    io = open(`bgzip -d -c $path`)
    (all_strains, chrom_rank, info_headers) = parse_vcf_header(io)
    pf = PeekedFile(io, "", false)
    advance!(pf)
    (pf, all_strains, chrom_rank, info_headers)
end

"""
    parse_vcf_record(line, n_samples) -> VCFRecord
"""
function parse_vcf_record(line::String, n_samples::Int)::VCFRecord
    fields = split(line, '\t')
    chrom = String(fields[1])
    pos   = parse(Int, fields[2])
    ref   = String(fields[4])
    alts  = String[String(a) for a in split(fields[5], ',')]
    info  = String(fields[8])
    fmt   = String(fields[9])
    format_keys = String[String(k) for k in split(fmt, ':')]
    sample_data = String[String(fields[9+i]) for i in 1:n_samples if 9+i <= length(fields)]
    VCFRecord(chrom, pos, ref, alts, info, format_keys, sample_data)
end

# Returns true if all ALTs have the same length as REF (SNP or MNP, not indel)
is_snp_record(record::VCFRecord) = all(length(a) == length(record.ref) for a in record.alts)

# Returns true if any ALT has the same length as REF. mergeVariantsByLocation.py's
# _recombine can merge decomposed rows back into one multi-ALT record that mixes a
# SNP allele with an indel allele (e.g. REF=A ALT=G,AT) — such a record is not a
# "SNP record" under is_snp_record, but still carries a real SNP allele.
has_snp_allele(record::VCFRecord) = any(length(a) == length(record.ref) for a in record.alts)

"""
    parse_format_field(format_keys, sample_str) -> Dict{String,String}
"""
function parse_format_field(format_keys::Vector{String}, sample_str::String)::Dict{String,String}
    result = Dict{String,String}()
    values = split(sample_str, ':')
    for (i, key) in enumerate(format_keys)
        i > length(values) && break
        result[key] = String(values[i])
    end
    result
end

# ---------------------------------------------------------------------------
# Coverage I/O (coverage.tsv produced by mergeCoverageBeds)
# ---------------------------------------------------------------------------

mutable struct CoverageFileHandle
    fh::IO
    sample_cols::Vector{Tuple{Int, String}}   # (col_index_1based, sample_name)
    peeked::Union{String, Nothing}
    exhausted::Bool
end

"""
    open_coverage_file(path) -> CoverageFileHandle

Opens coverage.tsv, reads the header to build a sample→column-index mapping,
and buffers the first data line.
"""
function open_coverage_file(path::String)::CoverageFileHandle
    fh = open(path, "r")
    header = readline(fh)
    fields = split(header, '\t')
    # Columns 4+ are sample names (1-based indexing: column 4 = index 4)
    sample_cols = Tuple{Int, String}[(i, String(fields[i])) for i in 4:length(fields)]
    first_line  = eof(fh) ? nothing : readline(fh)
    CoverageFileHandle(fh, sample_cols, first_line, first_line === nothing)
end

"""
    load_chrom_coverage!(cfh, chrom, chrom_rank, chrom_coverage)

Advances cfh past any lines for chromosomes that sort before `chrom` (by chrom_rank),
then reads all intervals for `chrom` into `chrom_coverage`, replacing prior contents.
Interval vectors are already sorted since coverage.tsv is position-sorted.
"""
function load_chrom_coverage!(
    cfh::CoverageFileHandle,
    chrom::String,
    chrom_rank::Dict{String, Int},
    chrom_coverage::Dict{String, Vector{Tuple{Int, Int, Float64}}}
)
    empty!(chrom_coverage)
    cfh.exhausted && return

    target_rank = get(chrom_rank, chrom, typemax(Int))

    # Advance past chromosomes that sort before the target
    while !cfh.exhausted
        fields    = split(cfh.peeked, '\t')
        line_rank = get(chrom_rank, String(fields[1]), typemax(Int))
        line_rank >= target_rank && break
        cfh.peeked    = eof(cfh.fh) ? nothing : readline(cfh.fh)
        cfh.exhausted = cfh.peeked === nothing
    end

    # Read all lines for this chrom
    while !cfh.exhausted
        fields = split(cfh.peeked, '\t')
        String(fields[1]) != chrom && break

        start_pos = parse(Int, fields[2])
        end_pos   = parse(Int, fields[3])

        for (col_idx, sample) in cfh.sample_cols
            col_idx > length(fields) && continue
            dp = parse(Float64, String(fields[col_idx]))
            if dp > 0.0
                push!(get!(chrom_coverage, sample, Tuple{Int,Int,Float64}[]),
                      (start_pos, end_pos, dp))
            end
        end

        cfh.peeked    = eof(cfh.fh) ? nothing : readline(cfh.fh)
        cfh.exhausted = cfh.peeked === nothing
    end
end

"""
    get_coverage(chrom_coverage, sample, pos) -> (covered, mean_dp)

Binary search for coverage at 0-based `pos`. Returns (false, 0.0) if not covered.
coverage.tsv uses 0-based half-open intervals [start, end) matching BED convention.
Pass VCF positions as `record.pos - 1` to convert from 1-based to 0-based.
"""
function get_coverage(
    chrom_coverage::Dict{String, Vector{Tuple{Int, Int, Float64}}},
    sample::String,
    pos::Int
)::Tuple{Bool, Float64}
    intervals = get(chrom_coverage, sample, nothing)
    (isnothing(intervals) || isempty(intervals)) && return (false, 0.0)
    idx = searchsortedlast(intervals, pos, by = x -> x[1])
    idx == 0 && return (false, 0.0)
    (_, end_, dp) = intervals[idx]
    return pos < end_ ? (true, dp) : (false, 0.0)
end

# ---------------------------------------------------------------------------
# TSV cache I/O
# Cache format: chrom \t pos \t transcript_id \t pos_in_cds  (one row per CDS transcript)
# ---------------------------------------------------------------------------

"""
    open_cache_peeked(path) -> PeekedFile

Opens the TSV cache file for reading. Returns a pre-exhausted PeekedFile
if the file is absent or empty. Skips comment lines (#-prefixed).
"""
function open_cache_peeked(path::String)::PeekedFile
    if !isfile(path) || filesize(path) == 0
        return PeekedFile(IOBuffer(""), "", true)
    end

    fh = open(path, "r")
    while !eof(fh)
        line = readline(fh)
        if !startswith(line, '#')
            pf = PeekedFile(fh, line, false)
            if isempty(line) && eof(fh)
                pf.exhausted = true
            end
            return pf
        end
    end
    close(fh)
    PeekedFile(IOBuffer(""), "", true)
end

"""
    parse_cache_tsv_record(line) -> (chrom, pos, transcript_id, pos_in_cds) or nothing
"""
function parse_cache_tsv_record(line::AbstractString)
    startswith(line, '#') && return nothing
    isempty(line)          && return nothing
    fields = split(line, '\t')
    length(fields) < 4 && return nothing
    chrom = String(fields[1])
    pos   = parse(Int, fields[2])
    tid   = String(fields[3])
    pic   = parse(Int, fields[4])
    (chrom, pos, tid, pic)
end


"""
    build_annotations_from_cache(entries, reference_strain, transcript_db, transcript_cache)
        -> Vector{PositionAnnotation}

Reconstructs PositionAnnotations from cached (transcript_id, pos_in_cds) pairs,
loading sequences from transcript_db as needed.
"""
function build_annotations_from_cache(
    entries::Vector{Tuple{String,Int}},
    reference_strain::String,
    transcript_db::SQLite.DB,
    transcript_cache::TranscriptSequenceCache
)::Vector{PositionAnnotation}
    isempty(entries) && return [PositionAnnotation(0, "", 0, 0, 0, "", "")]
    annotations = PositionAnnotation[]
    for (tid, pos_in_cds) in entries
        load_transcript_if_needed!(transcript_cache, transcript_db, tid)
        ref_seq = get_strain_seq(transcript_cache, tid, reference_strain)
        if isempty(ref_seq)
            ref_seq = get_strain_seq(transcript_cache, tid, "reference")
        end
        pos_in_codon_val = position_in_codon(pos_in_cds)
        ref_codon   = isempty(ref_seq) ? "" : extract_codon(ref_seq, pos_in_cds)
        ref_product = isempty(ref_codon) ? "" : translate_codon(ref_codon)
        push!(annotations, PositionAnnotation(1, tid, 0, pos_in_cds, pos_in_codon_val, ref_codon, ref_product))
    end
    annotations
end

function write_vcf_cache_header(fh::IO, all_strains::Vector{String}, info_headers::Vector{String})
    write(fh, "##fileformat=VCFv4.2\n")
    for h in info_headers
        write(fh, h, "\n")
    end
    write(fh, "##INFO=<ID=CANN,Number=.,Type=String,Description=\"Coding annotation entries, comma-separated. r-prefixed keys (r0,r1,...) = reference allele per transcript; k-prefixed keys (k0,k1,...) = alt allele per transcript. Format per entry: key|codon|aa|effect|transcript_id|pos_in_cds|pos_in_codon. Compound effects use '&' separator (e.g. missense&frameshift).\">\n")
    write(fh, "##FORMAT=<ID=CA,Number=1,Type=String,Description=\"CANN key(s) per GT allele. Alleles separated by '/' (unphased) or '|' (phased). Multiple transcript keys for one allele separated by ';'. 'r'=ref allele no CDS annotation, '.'=missing/no-call\">\n")
    write(fh, "##FORMAT=<ID=DFS,Number=1,Type=Integer,Description=\"Downstream of frameshift: 1 if this sample carries an upstream indel that disrupts the reading frame at this position, 0 otherwise.\">\n")
    chrom_line = join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", all_strains...], '\t')
    write(fh, chrom_line, "\n")
end

"""
    remap_sample_for_split(sample_str, format_keys, n_orig_alts, target_alt_i) -> String

When a multi-allelic record is split into one record per ALT, remap each sample's FORMAT
string so that GT allele indices are valid for a 1-ALT record:
  - target alt index (target_alt_i, 1-based in original) → 1
  - ref (0) → 0
  - any other alt index → 0 (treated as ref in this split record)
Also replaces GL with "." because GL has n*(n+1)/2 values for n alleles and SnpEff will
reject split records whose GL length no longer matches the single-ALT allele count.
"""
function remap_sample_for_split(sample_str::String, format_keys::Vector{String},
                                 n_orig_alts::Int, target_alt_i::Int)::String
    n_orig_alts == 1 && return sample_str   # nothing to remap for single-alt records
    values = split(sample_str, ':')
    result = String[String(v) for v in values]
    for (fi, key) in enumerate(format_keys)
        fi > length(result) && break
        if key == "GT"
            gt = result[fi]
            (isempty(gt) || gt == "." || gt == "./." || gt == ".|.") && continue
            sep_idx = findfirst(c -> c == '/' || c == '|', gt)
            sep = isnothing(sep_idx) ? nothing : gt[sep_idx]
            remap = idx -> idx == 0 ? 0 : (idx == target_alt_i ? 1 : 0)
            if isnothing(sep_idx)
                idx = tryparse(Int, gt)
                isnothing(idx) && continue
                result[fi] = string(remap(idx))
            else
                a1 = tryparse(Int, gt[1:sep_idx-1])
                a2 = tryparse(Int, gt[sep_idx+1:end])
                (isnothing(a1) || isnothing(a2)) && continue
                result[fi] = "$(remap(a1))$(sep)$(remap(a2))"
            end
        elseif key == "GL"
            result[fi] = "."
        end
    end
    join(result, ':')
end

function write_vcf_entry(fh::IO, chrom::String, pos::Int, ref::String, alt::String,
                          cann_str::String, gvcf_info::String,
                          format_keys::Vector{String}, sample_data::Vector{String},
                          ca_values::Vector{String}, dfs_values::Vector{String},
                          n_orig_alts::Int=1, target_alt_i::Int=1)
    info    = "$(gvcf_info);CANN=$(cann_str)"
    format  = join([format_keys..., "CA", "DFS"], ':')
    remapped = [remap_sample_for_split(s, format_keys, n_orig_alts, target_alt_i) for s in sample_data]
    samples = [i <= length(ca_values) ? "$(remapped[i]):$(ca_values[i]):$(get(dfs_values, i, "."))" : "$(remapped[i]):.:."
               for i in 1:length(remapped)]
    write(fh, join([chrom, string(pos), ".", ref, alt, ".", ".", info, format, samples...], '\t'), "\n")
end

# ---------------------------------------------------------------------------
# Variation construction from GVCF
# ---------------------------------------------------------------------------

"""
    gt_to_base(gt, ref, alts) -> String

Converts a GT field value to a base/allele string.
Returns "" for missing genotypes.
For diploid hets with SNP alleles, returns IUPAC ambiguity code.
"""
function gt_to_base(gt::String, ref::String, alts::Vector{String})::String
    (isempty(gt) || gt == "." || gt == "./." || gt == ".|.") && return ""

    sep_idx = findfirst(c -> c == '/' || c == '|', gt)

    if isnothing(sep_idx)
        # Haploid
        idx = parse(Int, gt)
        return idx == 0 ? ref : alts[idx]
    else
        a1_str = gt[1:sep_idx-1]
        a2_str = gt[sep_idx+1:end]
        (a1_str == "." || a2_str == ".") && return ""

        a1 = parse(Int, a1_str)
        a2 = parse(Int, a2_str)

        b1 = a1 == 0 ? ref : alts[a1]
        b2 = a2 == 0 ? ref : alts[a2]

        b1 == b2 && return b1  # homozygous

        # Het: IUPAC for single-char SNPs
        if length(b1) == 1 && length(b2) == 1
            iupac = get(IUPAC_COMPRESS, Set([b1[1], b2[1]]), nothing)
            !isnothing(iupac) && return string(iupac)
        end

        # Complex het: return the non-ref allele
        return a1 != 0 ? b1 : b2
    end
end

"""
    nonref_alt_alleles(gt, alts) -> Vector{String}

Returns the unique non-ref allele strings carried by this sample, in index order.
For GT="0/1" returns [alts[1]]; for GT="1/1" returns [alts[1]]; for GT="1/2"
returns [alts[1], alts[2]].  Returns [] for missing or ref-only GTs.
"""
function nonref_alt_alleles(gt::String, alts::Vector{String})::Vector{String}
    (isempty(gt) || gt == "." || gt == "./." || gt == ".|.") && return String[]
    sep_idx = findfirst(c -> c == '/' || c == '|', gt)
    idxs = if isnothing(sep_idx)
        [parse(Int, gt)]
    else
        a1 = gt[1:sep_idx-1]; a2 = gt[sep_idx+1:end]
        (a1 == "." || a2 == ".") && return String[]
        [parse(Int, a1), parse(Int, a2)]
    end
    seen = Set{Int}()
    result = String[]
    for i in idxs
        i == 0 && continue
        i in seen && continue
        push!(seen, i)
        i <= length(alts) && push!(result, alts[i])
    end
    result
end

"""
    compute_percent(fmt, gt) -> String

Computes AO/(RO+AO)*100 for the alt allele. For het calls (0/1), picks the non-ref allele
index so percent stores the alt fraction; ref fraction = 100 - percent.
Returns "0.0" for ref-only or missing data.
"""
function compute_percent(fmt::Dict{String,String}, gt::String)::String
    ao_str = get(fmt, "AO", "")
    ro_str = get(fmt, "RO", "0")
    isempty(ao_str) && return "0.0"

    sep_idx = findfirst(c -> c == '/' || c == '|', gt)
    if isnothing(sep_idx)
        aidx = parse(Int, gt)
        aidx == 0 && return "0.0"
    else
        a1 = parse(Int, gt[1:sep_idx-1])
        a2 = parse(Int, gt[sep_idx+1:end])
        aidx = a1 != 0 ? a1 : a2
        aidx == 0 && return "0.0"
    end

    ao_values = split(ao_str, ',')
    aidx > length(ao_values) && return "0.0"

    ao = parse(Float64, ao_values[aidx])
    ro = parse(Float64, isempty(ro_str) ? "0" : ro_str)
    total = ro + ao
    total <= 0 && return "0.0"

    @sprintf("%.2f", ao / total * 100)
end

"""
    build_variations_from_record(record, all_strains, undone_strains, chrom_coverage)
        -> Vector{Variation}

Builds per-strain Variation records from a VCF variant record.
For missing GTs, synthesizes a reference call if coverage.tsv shows the position covered.
"""
function build_variations_from_record(
    record::VCFRecord,
    all_strains::Vector{String},
    undone_strains::Set{String},
    chrom_coverage::Dict{String, Vector{Tuple{Int, Int, Float64}}},
    ploidy::Int=1
)::Vector{Variation}
    variations = Variation[]

    for (i, strain) in enumerate(all_strains)
        strain in undone_strains && continue
        i > length(record.sample_data) && continue

        fmt = parse_format_field(record.format_keys, record.sample_data[i])

        gt = get(fmt, "GT", "")
        if isempty(gt) || gt == "." || gt == "./." || gt == ".|."
            # No call: synthesize a reference Variation if position is covered
            (covered, dp) = get_coverage(chrom_coverage, strain, record.pos - 1)
            if covered
                v = Variation()
                v.sequence_source_id = record.chrom
                v.location           = record.pos
                v.strain             = strain
                v.reference          = record.ref
                v.base               = record.ref
                v.coverage           = string(dp)
                v.percent            = "100"
                v.quality            = "."
                v.pvalue             = "."
                v.snp_source_id      = "NGS_SNP.$(record.chrom).$(record.pos)"
                v.matches_reference  = 1
                v.ploidy             = ploidy
                push!(variations, v)
            end
            continue
        end

        dp_str = get(fmt, "DP", "0")
        dp     = isempty(dp_str) || dp_str == "." ? 0 : parse(Int, dp_str)

        base = gt_to_base(gt, record.ref, record.alts)
        if isempty(base) || base == "*"
            base == "*" && @warn "Unexpected * allele in GT at $(record.ref):$(record.pos) — should have been removed by mergeVariantsByLocation.py"
            continue
        end

        pct  = compute_percent(fmt, gt)
        gq   = get(fmt, "GQ", "0")

        # Detect het call: GT has separator and both alleles differ (e.g. 0/1)
        sep_idx  = findfirst(c -> c == '/' || c == '|', gt)
        is_het   = !isnothing(sep_idx) && gt[1:sep_idx-1] != gt[sep_idx+1:end]
        alt_alts = is_het ? nonref_alt_alleles(gt, record.alts) : String[]
        gt_ploidy = 1 + count(c -> c == '/' || c == '|', gt)

        v = Variation()
        v.sequence_source_id = record.chrom
        v.location           = record.pos
        v.strain             = strain
        v.reference          = record.ref
        v.base               = base
        v.alt_allele         = isempty(alt_alts) ? "" : alt_alts[1]
        v.coverage           = string(dp)
        v.percent            = pct
        v.quality            = gq
        v.pvalue             = "."
        v.snp_source_id      = "NGS_SNP.$(record.chrom).$(record.pos)"
        v.matches_reference  = (base == record.ref) ? 1 : 0
        v.ploidy             = gt_ploidy

        push!(variations, v)
    end

    variations
end

# ---------------------------------------------------------------------------
# Sorted-merge helpers: sort keys over VCF lines
# ---------------------------------------------------------------------------

"""
    peek_sort_key(line, chrom_rank) -> (Int, Int)

Returns (chrom_rank, pos) from a VCF data line for sorted-merge ordering.
"""
function peek_sort_key(line::String, chrom_rank::Dict{String,Int})::Tuple{Int,Int}
    tab1 = findfirst('\t', line)
    isnothing(tab1) && return (typemax(Int), typemax(Int))
    chrom = line[1:tab1-1]
    rest  = line[tab1+1:end]
    tab2  = findfirst('\t', rest)
    pos_str = isnothing(tab2) ? rest : rest[1:tab2-1]
    rank = get(chrom_rank, chrom, typemax(Int))
    (rank, parse(Int, pos_str))
end


# ---------------------------------------------------------------------------
# Resource management
# ---------------------------------------------------------------------------

"""
    initialize_processing_context(args, all_strains) -> ProcessingContext
"""
function initialize_processing_context(args, all_strains::Vector{String})
    undone_strains = Set{String}()
    undone_strains_file = get(args, "undone_strains_file", "")
    if !isempty(undone_strains_file) && isfile(undone_strains_file)
        open(undone_strains_file, "r") do fh
            for line in eachline(fh)
                s = strip(line)
                !isempty(s) && push!(undone_strains, s)
            end
        end
    end

    (cds_intervals, transcript_info) = parse_gtf(args["gtf_file"])

    transcript_db = SQLite.DB(args["transcript_db"])
    indel_db      = SQLite.DB(args["indel_db"])

    fs_info = precompute_frameshifts(indel_db, transcript_info)

    ploidy = haskey(args, "ploidy") ? parse(Int, args["ploidy"]) : 1

    ProcessingContext(
        args["reference_strain"],
        undone_strains,
        cds_intervals,
        transcript_info,
        transcript_db,
        indel_db,
        fs_info,
        all_strains,
        Dict{String,Int}(name => i for (i, name) in enumerate([all_strains..., args["reference_strain"]])),
        ploidy
    )
end

"""
    open_output_writers(output_vcf, reference_strain, all_strains) -> OutputWriters
"""
function open_output_writers(output_vcf::String, reference_strain::String,
                              all_strains::Vector{String})
    vcf_fh = open(output_vcf, "w")
    snp_fh = open("snpFeature.dat", "w")
    write(snp_fh, "location\tseq_id\treference_strain\tref_allele\tmajor_allele\tminor_allele\tmajor_allele_strain_count\tminor_allele_strain_count\tmajor_allele_frequency\tminor_allele_frequency\tdistinct_strain_count\tdistinct_allele_count\ttotal_ploidy_count\tis_coding\tvariant_type\tmajor_differs_from_reference\tis_singleton\thet_strain_count\tcalled_strain_count\tno_call_strain_count\tcall_rate\n")
    allele_fh = open("allele.dat", "w")
    write(allele_fh, "location\tseq_id\tallele\tdistinct_strain_count\tallele_frequency\tavg_coverage\tavg_percent\tstrain_ids\tmatches_reference\n")
    tp_fh = open("transcript_product.dat", "w")
    write(tp_fh, "#seq_id\tlocation\ttranscript_id\tpos_in_cds\tpos_in_protein\tcodon\tpos_in_codon\tcount\tproduct\tmatches_ref_codon\tmatches_ref_product\tdownstream_of_frameshift_strain_ids\n")
    hsss = open_hsss_writers(reference_strain, all_strains)
    OutputWriters(vcf_fh, snp_fh, allele_fh, tp_fh, hsss)
end

function close_processing_context(ctx::ProcessingContext)
    close(ctx.transcript_db)
    close(ctx.indel_db)
end

function close_output_writers(writers::OutputWriters)
    close(writers.vcf_fh)
    close(writers.snp_fh)
    close(writers.allele_fh)
    close(writers.transcript_product_fh)
    close_hsss_writers(writers.hsss)
end

# ---------------------------------------------------------------------------
# Variation / position checking
# ---------------------------------------------------------------------------

function has_variation(variations::Vector{Variation})
    length(Set(v.base for v in variations)) > 1
end

# ---------------------------------------------------------------------------
# Annotation logic
# ---------------------------------------------------------------------------

"""
    annotate_position_all(seq_id, location, ctx, transcript_cache) -> Vector{PositionAnnotation}

Returns one PositionAnnotation per CDS transcript overlapping location.
Returns a single non-coding annotation if location is intergenic.
"""
function annotate_position_all(
    seq_id::String,
    location::Int,
    ctx::ProcessingContext,
    transcript_cache::TranscriptSequenceCache
)::Vector{PositionAnnotation}
    cds_hits = find_all_cds(ctx.cds_intervals, seq_id, location)

    if isempty(cds_hits)
        return [PositionAnnotation(0, "", 0, 0, 0, "", "")]
    end

    annotations = PositionAnnotation[]
    seen_transcripts = Set{String}()

    for cds_hit in cds_hits
        transcript_id = cds_hit.transcript_id
        transcript_id in seen_transcripts && continue
        push!(seen_transcripts, transcript_id)

        tinfo            = ctx.transcript_info[transcript_id]
        cds_number       = cds_hit.cds_number
        pos_in_cds       = compute_position_in_cds(tinfo, location, cds_number)
        pos_in_codon_val = position_in_codon(pos_in_cds)

        debug_log("    Loading sequences for transcript: ", transcript_id)
        load_transcript_if_needed!(transcript_cache, ctx.transcript_db, transcript_id)

        ref_seq = get_strain_seq(transcript_cache, transcript_id, ctx.reference_strain)
        if isempty(ref_seq)
            ref_seq = get_strain_seq(transcript_cache, transcript_id, "reference")
        end

        ref_codon   = ""
        ref_product = ""
        if !isempty(ref_seq)
            ref_codon   = extract_codon(ref_seq, pos_in_cds)
            ref_product = translate_codon(ref_codon)
        end

        push!(annotations, PositionAnnotation(1, transcript_id, cds_number, pos_in_cds,
                                              pos_in_codon_val, ref_codon, ref_product))
    end

    annotations
end

"""
    annotate_variations!(variations, annotation, ctx, transcript_cache)

Fills per-strain codon/product data into each Variation record.
"""
function annotate_variations!(
    variations::Vector{Variation},
    annotation::PositionAnnotation,
    ctx::ProcessingContext,
    transcript_cache::TranscriptSequenceCache
)
    for v in variations
        v.is_coding  = annotation.is_coding
        v.transcript = annotation.transcript_id
        v.cds_number = annotation.cds_number

        annotation.is_coding == 1 || continue

        v.position_in_cds    = annotation.pos_in_cds
        v.position_in_codon  = annotation.pos_in_codon_val
        v.reference_codon    = annotation.ref_codon

        strain_seq = get_strain_seq(transcript_cache, annotation.transcript_id, v.strain)
        if !isempty(strain_seq)
            is_fs = check_downstream_of_frameshift(ctx.fs_info, v.strain, annotation.transcript_id, annotation.pos_in_cds) == 1
            if is_fs
                v.codon   = "."
                v.product = String[]
            else
                shift        = get_indel_shift(ctx.indel_db, annotation.transcript_id, v.strain, annotation.pos_in_cds)
                adjusted_pos = annotation.pos_in_cds + shift
                strain_codon = extract_codon(strain_seq, adjusted_pos)
                v.codon   = strain_codon
                expanded  = expand_codon(strain_codon)
                v.product = [translate_codon(c) for c in expanded]
            end
        else
            v.codon   = "NNN"
            v.product = ["X"]
        end
    end
end

"""
    build_reference_variation(variations, annotation, seq_id, location, reference_strain) -> Variation
"""
function build_reference_variation(
    variations::Vector{Variation},
    annotation::PositionAnnotation,
    seq_id::String,
    location::Int,
    reference_strain::String
)
    ref_allele = variations[1].reference
    isempty(ref_allele) && (ref_allele = variations[1].base)

    adjacent_snp_causes_product_difference = 0
    if annotation.is_coding == 1
        for v in variations
            product_str = join(v.product, ":")
            if product_str != annotation.ref_product
                adjacent_snp_causes_product_difference = 1
                break
            end
        end
    end

    for v in variations
        if annotation.is_coding == 1
            v.has_nonsynonomous = (join(v.product, ":") != annotation.ref_product) ? 1 : 0
        end
    end

    Variation(
        seq_id, location, reference_strain,
        ref_allele, ref_allele,
        "",      # alt_allele: empty — reference variation is never a het call
        "", "",
        "", "",
        "",
        annotation.is_coding,
        annotation.pos_in_cds,
        annotation.pos_in_codon_val,
        0,
        annotation.transcript_id,
        annotation.is_coding == 1 ? [annotation.ref_product] : String[],
        annotation.ref_codon,
        annotation.ref_codon,
        adjacent_snp_causes_product_difference,
        annotation.cds_number,
        1,   # matches_reference
        1    # ploidy: reference genome is a single representative
    )
end

# ---------------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------------

"""
    compute_allele_weight_map(variations) -> (Dict{String,Int}, Int)

Returns (allele_weight_counts, total_weight).  Het calls contribute weight 1
to each of ref and alt; hom calls contribute v.ploidy to v.base.
"""
function compute_allele_weight_map(variations::Vector{Variation})::Tuple{Dict{String,Int}, Int}
    weights = Dict{String,Int}()
    total   = 0
    for v in variations
        if !isempty(v.alt_allele)
            weights[v.reference]  = get(weights, v.reference,  0) + 1
            weights[v.alt_allele] = get(weights, v.alt_allele, 0) + 1
            total += 2
        else
            weights[v.base] = get(weights, v.base, 0) + v.ploidy
            total += v.ploidy
        end
    end
    (weights, total)
end

function write_snp_feature(
    snp_fh::IO,
    variations::Vector{Variation},
    is_coding::Int,
    seq_id::String,
    location::Int,
    reference_strain::String,
    sequenced_strains::Vector{String}
)
    ref_allele = variations[1].reference
    isempty(ref_allele) && (ref_allele = variations[1].base)

    (allele_weights, total_ploidy_count) = compute_allele_weight_map(variations)

    # strain count is per-strain, not ploidy-weighted; compute separately
    allele_counts = Dict{String,Int}()
    strain_set    = Set{String}()
    for v in variations
        if !isempty(v.alt_allele)
            allele_counts[v.reference]  = get(allele_counts, v.reference,  0) + 1
            allele_counts[v.alt_allele] = get(allele_counts, v.alt_allele, 0) + 1
        else
            allele_counts[v.base] = get(allele_counts, v.base, 0) + 1
        end
        push!(strain_set, v.strain)
    end

    distinct_strain_count = length(strain_set)
    distinct_allele_count = length(allele_counts)

    n_alt_alleles  = count(a -> a != ref_allele, keys(allele_counts))
    sorted_alleles = sort(collect(keys(allele_counts));
                         by = a -> (n_alt_alleles >= 2 && a == ref_allele ? 1 : 0,
                                    -allele_counts[a], a))

    major_allele              = sorted_alleles[1]
    minor_allele              = length(sorted_alleles) > 1 ? sorted_alleles[2] : ""
    major_allele_strain_count = allele_counts[major_allele]
    minor_allele_strain_count = length(sorted_alleles) > 1 ? allele_counts[minor_allele] : ""
    major_allele_frequency    = @sprintf("%.4f", allele_weights[major_allele] / total_ploidy_count)
    minor_allele_frequency    = length(sorted_alleles) > 1 ?
        @sprintf("%.4f", allele_weights[minor_allele] / total_ploidy_count) : ""

    # --- precompute columns -------------------------------------------------
    # variant_type: classify the alt alleles present at this locus by length
    # relative to the reference. A same-length alt is a substitution (SNV);
    # a different-length alt is an indel. Both present => MIXED.
    has_snp   = false
    has_indel = false
    for v in variations
        if !isempty(v.alt_allele)                       # het call: alt vs ref
            length(v.reference) == length(v.alt_allele) ? (has_snp = true) : (has_indel = true)
        elseif v.base != v.reference                    # hom alt call: base vs ref
            length(v.reference) == length(v.base) ? (has_snp = true) : (has_indel = true)
        end
    end
    variant_type = has_snp && has_indel ? "MIXED" : (has_indel ? "INDEL" : "SNV")

    major_differs_from_reference = major_allele != ref_allele ? "1" : "0"
    is_singleton = (minor_allele_strain_count isa Integer && minor_allele_strain_count == 1) ? "1" : "0"
    het_strain_count = count(v -> !isempty(v.alt_allele), variations)

    # call rate is over sequenced samples only; the reference is a synthetic
    # always-present strain and is excluded from both numerator and denominator.
    called_strain_count  = count(s -> s != reference_strain, strain_set)
    total_sequenced      = length(sequenced_strains)
    no_call_strain_count = max(0, total_sequenced - called_strain_count)
    call_rate = total_sequenced > 0 ?
        @sprintf("%.4f", called_strain_count / total_sequenced) : ""

    write(snp_fh, join([
        string(location),
        seq_id,
        reference_strain,
        ref_allele,
        major_allele,
        minor_allele,
        string(major_allele_strain_count),
        string(minor_allele_strain_count),
        major_allele_frequency,
        minor_allele_frequency,
        string(distinct_strain_count),
        string(distinct_allele_count),
        string(total_ploidy_count),
        string(is_coding),
        variant_type,
        major_differs_from_reference,
        is_singleton,
        string(het_strain_count),
        string(called_strain_count),
        string(no_call_strain_count),
        call_rate
    ], "\t"), "\n")
end

function write_allele_file(
    allele_fh::IO,
    variations::Vector{Variation},
    seq_id::String,
    location::Int,
    sample_id_map::Dict{String,Int}
)
    # allele_entries: allele -> [(strain, coverage, percent)]
    # allele_weights comes from compute_allele_weight_map for frequency denominator
    allele_entries = Dict{String, Vector{Tuple{String, Float64, Float64}}}()

    for v in variations
        cov = isempty(v.coverage) ? 0.0 : parse(Float64, v.coverage)
        pct = isempty(v.percent)  ? 0.0 : parse(Float64, v.percent)

        if !isempty(v.alt_allele)
            push!(get!(allele_entries, v.reference,  []), (v.strain, cov, 100.0 - pct))
            push!(get!(allele_entries, v.alt_allele, []), (v.strain, cov, pct))
        else
            push!(get!(allele_entries, v.base, []), (v.strain, cov, pct))
        end
    end

    (allele_weights, total_weight) = compute_allele_weight_map(variations)
    ref_allele = variations[1].reference

    for (allele, entries) in allele_entries
        distinct_strains = Set{String}()
        strain_ids       = Set{Int}()
        sum_coverage     = 0.0
        sum_percent      = 0.0

        for (strain, cov, pct) in entries
            push!(distinct_strains, strain)
            sid = get(sample_id_map, strain, 0)
            if sid > 0
                push!(strain_ids, sid)
            else
                @warn "strain not found in sample_id_map, omitted from strain_ids" strain
            end
            sum_coverage += cov
            sum_percent  += pct
        end

        n       = length(entries)
        ids_str = "{" * join(sort(collect(strain_ids)), ",") * "}"
        matches_ref = allele == ref_allele ? 1 : 0
        write(allele_fh, join([
            string(location),
            seq_id,
            allele,
            string(length(distinct_strains)),
            @sprintf("%.4f", allele_weights[allele] / total_weight),
            @sprintf("%.2f", sum_coverage / n),
            @sprintf("%.2f", sum_percent  / n),
            ids_str,
            string(matches_ref)
        ], "\t"), "\n")
    end

end

function write_transcript_product(
    fh::IO,
    variations::Vector{Variation},
    annotation::PositionAnnotation,
    seq_id::String,
    location::Int,
    sample_id_map::Dict{String,Int}
)
    annotation.is_coding != 1 && return

    pos_in_protein = div(annotation.pos_in_cds - 1, 3) + 1

    dfs_ids = sort([sample_id_map[v.strain] for v in variations
                    if v.downstream_of_frameshift == 1 && haskey(sample_id_map, v.strain)])
    dfs_str = isempty(dfs_ids) ? "" : "{" * join(dfs_ids, ",") * "}"

    all_product_counts = Dict{String,Int}()
    for v in variations
        for p in v.product
            all_product_counts[p] = get(all_product_counts, p, 0) + 1
        end
    end

    seen_codons = Set{String}()
    for v in variations
        isempty(v.codon) && continue
        v.downstream_of_frameshift == 1 && continue
        for ec in expand_codon(v.codon)
            push!(seen_codons, ec)
        end
    end

    for ec in seen_codons
        product             = translate_codon(ec)
        count               = get(all_product_counts, product, 0)
        matches_ref_codon   = ec == annotation.ref_codon ? 1 : 0
        matches_ref_product = product == annotation.ref_product ? 1 : 0
        write(fh, join([
            seq_id,
            string(location),
            annotation.transcript_id,
            string(annotation.pos_in_cds),
            string(pos_in_protein),
            ec,
            string(annotation.pos_in_codon_val),
            string(count),
            product,
            string(matches_ref_codon),
            string(matches_ref_product),
            dfs_str
        ], "\t"), "\n")
    end
end

# ---------------------------------------------------------------------------
# CA FORMAT field helpers
# ---------------------------------------------------------------------------

"""
    gt_to_ca(gt, this_alt_idx, this_alt_keys, ref_keys) -> String

Maps a sample's GT string to a CA allele-key string.
- `this_alt_keys`: CANN k-keys for this cache line's alt (one per transcript); ';'-joined.
- `ref_keys`: CANN r-keys for the reference allele (one per transcript); ';'-joined.
  When ref_keys is empty (non-coding position) 'r' is emitted as a plain marker.
Preserves the phasing separator ('/' or '|') from the original GT.
'.' is emitted only for missing/no-call genotypes.
"""
function gt_to_ca(gt::String, this_alt_idx::Union{Int,Nothing},
                  this_alt_keys::Vector{String}, ref_keys::Vector{String})::String
    (isempty(gt) || gt == "." || gt == "./." || gt == ".|.") && return "."

    alt_key_str = isempty(this_alt_keys) ? "." : join(this_alt_keys, ';')
    ref_key_str = isempty(ref_keys)      ? "r" : join(ref_keys, ';')

    sep_idx = findfirst(c -> c == '/' || c == '|', gt)

    allele_to_ca(a_str) = begin
        idx = tryparse(Int, a_str)
        isnothing(idx) && return "."
        idx == 0 && return ref_key_str
        (this_alt_idx !== nothing && idx == this_alt_idx) ? alt_key_str : "."
    end

    if isnothing(sep_idx)
        return allele_to_ca(gt)
    else
        sep_char = gt[sep_idx]
        ca1 = allele_to_ca(gt[1:sep_idx-1])
        ca2 = allele_to_ca(gt[sep_idx+1:end])
        return "$(ca1)$(sep_char)$(ca2)"
    end
end

"""
    build_ca_values(format_keys, sample_data, all_alts, this_alt, this_alt_keys, ref_keys) -> Vector{String}

Returns one CA string per sample. `this_alt_keys` lists the k-keys for this alt's
transcripts; `ref_keys` lists the r-keys for the reference allele's transcripts.
Multiple keys for one allele are ';'-joined.
"""
function build_ca_values(
    format_keys::Vector{String},
    sample_data::Vector{String},
    all_alts::Vector{String},
    this_alt::String,
    ref_keys::Vector{String},
    all_strains::Vector{String},
    strain_to_alt_keys::Dict{String, Vector{String}}
)::Vector{String}
    this_alt_idx = findfirst(==(this_alt), all_alts)
    result = String[]
    for (i, sd) in enumerate(sample_data)
        strain = i <= length(all_strains) ? all_strains[i] : ""
        gt     = get(parse_format_field(format_keys, sd), "GT", ".")
        per_sample_keys = get(strain_to_alt_keys, strain, String[])
        push!(result, gt_to_ca(gt, this_alt_idx, per_sample_keys, ref_keys))
    end
    result
end

const AA_THREE_LETTER = Dict(
    "A"=>"Ala","R"=>"Arg","N"=>"Asn","D"=>"Asp","C"=>"Cys",
    "Q"=>"Gln","E"=>"Glu","G"=>"Gly","H"=>"His","I"=>"Ile",
    "L"=>"Leu","K"=>"Lys","M"=>"Met","F"=>"Phe","P"=>"Pro",
    "S"=>"Ser","T"=>"Thr","W"=>"Trp","Y"=>"Tyr","V"=>"Val",
    "*"=>"Ter",
)

"""
    substitution_hgvs(pos_in_cds, ref_codon, alt_codon, pic, ref_aa, alt_aa) -> (String, String)

Builds (HGVS.c, HGVS.p) for a single-base coding substitution. Bases are read
from the strand-oriented codons at position `pic`, so no separate strand
handling is required. Returns "." for either string it cannot form (ambiguous
or non-triplet codon, unknown amino acid, or no base change at `pic`).
"""
function substitution_hgvs(pos_in_cds::Int, ref_codon::String, alt_codon::String,
                           pic::Int, ref_aa::String, alt_aa::String)::Tuple{String,String}
    hgvs_c = "."
    if length(ref_codon) == 3 && length(alt_codon) == 3 &&
       !occursin(r"[^ACGTacgt]", ref_codon) && !occursin(r"[^ACGTacgt]", alt_codon) &&
       1 <= pic <= 3
        rb = ref_codon[pic]
        ab = alt_codon[pic]
        if rb != ab
            hgvs_c = "c.$(pos_in_cds)$(rb)>$(ab)"
        end
    end

    hgvs_p = "."
    ref3 = get(AA_THREE_LETTER, ref_aa, "")
    alt3 = get(AA_THREE_LETTER, alt_aa, "")
    if !isempty(ref3) && !isempty(alt3)
        protpos = div(pos_in_cds - 1, 3) + 1
        hgvs_p = if ref_aa == alt_aa
            "p.$(ref3)$(protpos)="
        elseif protpos == 1 && ref_aa == "M"
            "p.Met1?"
        elseif alt_aa == "*"
            "p.$(ref3)$(protpos)Ter"
        else
            "p.$(ref3)$(protpos)$(alt3)"
        end
    end

    (hgvs_c, hgvs_p)
end

# ---------------------------------------------------------------------------
# CANN annotation
# ---------------------------------------------------------------------------

"""
    build_ref_cann_entry(key, annotation) -> String

Builds the r-keyed CANN entry for the reference allele at a coding position.
Format mirrors alt entries: key|codon|aa|effect|transcript_id|pos_in_cds|pos_in_codon
"""
function build_ref_cann_entry(key::String, annotation::PositionAnnotation)::String
    annotation.is_coding != 1 && return "."
    codon = isempty(annotation.ref_codon)   ? "." : annotation.ref_codon
    aa    = isempty(annotation.ref_product) ? "." : annotation.ref_product
    "$(key)|$(codon)|$(aa)|reference|$(annotation.transcript_id)|$(annotation.pos_in_cds)|$(annotation.pos_in_codon_val)"
end

"""
    build_cann_string(ref_allele, alt_allele, v, annotation) -> String

Builds the CANN entry for one (ref, alt) pair.
Returns "." for non-coding positions.
"""
function build_cann_string(
    ref_allele::String,
    alt_allele::String,
    v::Variation,
    annotation::PositionAnnotation
)::String
    annotation.is_coding != 1 && return "."

    tid        = annotation.transcript_id
    pos_in_cds = annotation.pos_in_cds
    pic        = annotation.pos_in_codon_val

    ref_len = length(ref_allele)
    alt_len = length(alt_allele)
    is_indel = (ref_len != alt_len)

    # A complex variant has both a substitution and an indel component.
    # FreeBayes emits: [shared prefix][variant bases][inserted/deleted bases] — no shared suffix.
    # After stripping the common prefix, if both tails are non-empty the variant is complex.
    # This generalises the old anchor-only check (ref[1] != alt[1]) to handle cases like
    # ATA→ACCTG where the anchor matches but downstream bases differ before the net insert.
    is_complex = if is_indel && !isempty(ref_allele) && !isempty(alt_allele)
        pfx = 0
        for i in 1:min(ref_len, alt_len)
            ref_allele[i] == alt_allele[i] || break
            pfx += 1
        end
        !isempty(ref_allele[pfx+1:end]) && !isempty(alt_allele[pfx+1:end])
    else
        false
    end

    is_pure_indel = is_indel && !is_complex

    if is_pure_indel
        len_diff = alt_len - ref_len
        structural = if len_diff % 3 != 0
            "frameshift"
        elseif len_diff > 0
            "inframe_insertion"
        else
            "inframe_deletion"
        end
        return "k0|.|.|$(structural)|$(tid)|$(pos_in_cds)|$(pic)"
    end

    # SNP or complex variant: compute amino acid effect
    codon       = isempty(v.codon) ? "." : v.codon
    unique_prods = unique(v.product)
    product_str  = isempty(unique_prods) ? "." : join(unique_prods, "/")

    # Codon/product suppressed because strain is downstream of a frameshift
    if codon == "." && isempty(unique_prods)
        return "k0|.|.|downstream_frameshift|$(tid)|$(pos_in_cds)|$(pic)"
    end

    # Codon contains ambiguous base(s) — skip product and effect
    if occursin(r"[NnXx]", codon)
        return "k0|$(codon)|.|.|$(tid)|$(pos_in_cds)|$(pic)"
    end

    has_stop = any(p == "*" for p in unique_prods)
    aa_effect = if has_stop
        "nonsense"
    elseif product_str != annotation.ref_product
        "missense"
    else
        "synonymous"
    end

    if !is_indel
        return "k0|$(codon)|$(product_str)|$(aa_effect)|$(tid)|$(pos_in_cds)|$(pic)"
    else
        # Complex: indel with SNP at anchor position
        len_diff = alt_len - ref_len
        structural = if len_diff % 3 != 0
            "frameshift"
        elseif len_diff > 0
            "inframe_insertion"
        else
            "inframe_deletion"
        end
        return "k0|$(codon)|$(product_str)|$(aa_effect)&$(structural)|$(tid)|$(pos_in_cds)|$(pic)"
    end
end


# ---------------------------------------------------------------------------
# Core record handler — extracted helpers
# ---------------------------------------------------------------------------

"""
    pick_snp_record(records) -> VCFRecord

Returns the first SNP record from a group of records sharing the same position.
Falls back to records[1] when no SNP record exists (all-indel group).
"""
function pick_snp_record(records::Vector{VCFRecord})::VCFRecord
    snp_records = filter(is_snp_record, records)
    isempty(snp_records) ? records[1] : snp_records[1]
end

"""
    build_ref_cann_entries(annotations) -> (Vector{String}, Vector{String})

Builds the r-keyed CANN strings for the reference allele.
Returns (ref_keys, ref_cann_entries). Non-coding and dot entries are skipped.
"""
function build_ref_cann_entries(annotations::Vector{PositionAnnotation})::Tuple{Vector{String}, Vector{String}}
    ref_keys         = String[]
    ref_cann_entries = String[]
    for annotation in annotations
        annotation.is_coding == 0 && continue
        key   = "r$(length(ref_keys))"
        entry = build_ref_cann_entry(key, annotation)
        entry == "." && continue
        push!(ref_keys, key)
        push!(ref_cann_entries, entry)
    end
    (ref_keys, ref_cann_entries)
end

"""
    fill_missing_coverage_gt(record, all_strains, chrom_coverage) -> Vector{String}

Returns a copy of record.sample_data with GT=0 and DP filled in for samples
that have a missing/dot GT but are covered at this position.
"""
function fill_missing_coverage_gt(
    record::VCFRecord,
    all_strains::Vector{String},
    chrom_coverage::Dict{String, Vector{Tuple{Int, Int, Float64}}},
    ploidy::Int=1
)::Vector{String}
    gt_idx  = findfirst(==("GT"), record.format_keys)
    dp_idx  = findfirst(==("DP"), record.format_keys)
    gt_fill = join(fill("0", ploidy), "/")
    modified = copy(record.sample_data)
    for (i, strain) in enumerate(all_strains)
        i > length(modified) && continue
        fmt = parse_format_field(record.format_keys, modified[i])
        gt = get(fmt, "GT", "")
        (isempty(gt) || gt == "." || gt == "./." || gt == ".|.") || continue
        (covered, dp) = get_coverage(chrom_coverage, strain, record.pos - 1)
        covered || continue
        fields = fill(".", length(record.format_keys))
        !isnothing(gt_idx) && (fields[gt_idx] = gt_fill)
        !isnothing(dp_idx) && (fields[dp_idx] = string(round(Int, dp)))
        modified[i] = join(fields, ":")
    end
    modified
end

"""
    assign_cann_keys(alt_strain_entries, all_strains) -> (Dict, Dict)

Deduplicates CANN strings and assigns k0/k1/... keys in canonical strain order.
Returns (alt_cann_entries, alt_strain_to_ca):
  - alt_cann_entries:  alt -> ordered unique keyed CANN strings
  - alt_strain_to_ca:  alt -> strain -> assigned key vector
"""
function assign_cann_keys(
    alt_strain_entries::Dict{String, Dict{String, Vector{String}}},
    all_strains::Vector{String}
)::Tuple{Dict{String, Vector{String}}, Dict{String, Dict{String, Vector{String}}}}
    alt_cann_entries = Dict{String, Vector{String}}()
    alt_strain_to_ca = Dict{String, Dict{String, Vector{String}}}()

    for (alt, strain_map) in alt_strain_entries
        entry_to_key  = Dict{String, String}()
        ordered_keyed = String[]

        for strain in all_strains
            strain_entries = get(strain_map, strain, nothing)
            isnothing(strain_entries) && continue
            for entry in strain_entries
                entry == "." && continue
                if !haskey(entry_to_key, entry)
                    key = "k$(length(entry_to_key))"
                    entry_to_key[entry] = key
                    push!(ordered_keyed, replace(entry, r"^k0\|" => "$(key)|", count=1))
                end
            end
        end

        alt_cann_entries[alt] = ordered_keyed

        strain_keys = Dict{String, Vector{String}}()
        for (strain, strain_entries) in strain_map
            strain_keys[strain] = [get(entry_to_key, e, ".") for e in strain_entries]
        end
        alt_strain_to_ca[alt] = strain_keys
    end

    (alt_cann_entries, alt_strain_to_ca)
end

"""
    collect_cann_entries_for_annotation(variations, annotation, record, reference_strain, all_strains, strain_idx_map)
        -> Dict{String, Dict{String, Vector{String}}}

Collects per-alt, per-strain CANN entry strings for one annotation (transcript).
Skips the reference strain and variations that match the reference.
Returns alt -> strain -> [entry, ...].
"""
function collect_cann_entries_for_annotation(
    variations::Vector{Variation},
    annotation::PositionAnnotation,
    record::VCFRecord,
    reference_strain::String,
    all_strains::Vector{String},
    strain_idx_map::Dict{String, Int}
)::Dict{String, Dict{String, Vector{String}}}
    result = Dict{String, Dict{String, Vector{String}}}()
    for v in variations
        v.strain == reference_strain && continue
        v.matches_reference == 1 && continue
        sidx = get(strain_idx_map, v.strain, 0)
        (sidx == 0 || sidx > length(record.sample_data)) && continue
        fmt = parse_format_field(record.format_keys, record.sample_data[sidx])
        gt  = get(fmt, "GT", "")
        for alt_allele in nonref_alt_alleles(gt, record.alts)
            if alt_allele == "*"
                @warn "Unexpected * allele in CANN annotation at $(record.ref):$(record.pos) — should have been removed by mergeVariantsByLocation.py"
                continue
            end
            entry = build_cann_string(record.ref, alt_allele, v, annotation)
            strain_map = get!(result, alt_allele, Dict{String, Vector{String}}())
            push!(get!(strain_map, v.strain, String[]), entry)
        end
    end
    result
end

# ---------------------------------------------------------------------------
# Core record handler
# ---------------------------------------------------------------------------

"""
    handle_variant_record!(records, cache_entries, ctx, writers, transcript_cache, all_strains, chrom_coverage) -> Bool

Processes one variant VCF record end-to-end. Returns true if output was written.
cache_entries: Vector of (transcript_id, pos_in_cds) pairs from the TSV cache, or empty.
"""
function handle_variant_record!(
    records::Vector{VCFRecord},
    cache_entries::Vector{Tuple{String,Int}},
    ctx::ProcessingContext,
    writers::OutputWriters,
    transcript_cache::TranscriptSequenceCache,
    all_strains::Vector{String},
    chrom_coverage::Dict{String, Vector{Tuple{Int, Int, Float64}}}
)::Bool
    seq_id   = records[1].chrom
    location = records[1].pos
    debug_log("Processing position: ", seq_id, ":", location)

    # Build variations from all records at this position (SNPs + indels combined)
    variations = Variation[]
    for r in records
        append!(variations, build_variations_from_record(r, all_strains, ctx.undone_strains, chrom_coverage, ctx.ploidy))
    end
    isempty(variations) && return false

    # For CANN annotation and VCF output: use only the first SNP record.
    # Indel records at the same position are included in variation tables but not AA annotation.
    record = pick_snp_record(records)

    # Determine annotations: one per overlapping transcript (cache or fresh GTF lookup)
    annotations = if !isempty(cache_entries)
        debug_log("    Cache hit at ", seq_id, ":", location)
        build_annotations_from_cache(cache_entries, ctx.reference_strain, ctx.transcript_db, transcript_cache)
    else
        annotate_position_all(seq_id, location, ctx, transcript_cache)
    end

    # Accumulate per-sample CANN entries across all annotations (transcripts).
    # alt_strain_entries: alt -> strain -> [entry per transcript, in annotation order]
    alt_strain_entries = Dict{String, Dict{String, Vector{String}}}()
    first_annotation   = annotations[1]
    any_output         = false
    first_all_vars        = nothing
    # Collect products from non-ref strain variations across all transcript annotations.
    # Reference strain product (annotation.ref_product) is intentionally excluded —
    # product_code reflects what the non-ref variants encode.
    all_annotation_products = String[]

    # Map strain name -> sample index for GT lookup when keying CANN entries by original alt allele
    strain_idx_map = Dict{String, Int}(s => i for (i, s) in enumerate(all_strains))

    # allele.dat is per-position (not per-transcript): written once after the first annotation pass
    allele_written = false

    for annotation in annotations
        annotate_variations!(variations, annotation, ctx, transcript_cache)

        for v in variations
            if !isempty(v.transcript)
                v.downstream_of_frameshift = check_downstream_of_frameshift(
                    ctx.fs_info, v.strain, v.transcript, v.position_in_cds)
            end
        end

        ref_variation = build_reference_variation(variations, annotation, seq_id, location, ctx.reference_strain)
        all_vars = vcat(variations, [ref_variation])

        has_variation(all_vars) || continue
        any_output = true

        if isnothing(first_all_vars)
            first_all_vars = all_vars
        end

        for v in variations
            append!(all_annotation_products, v.product)
        end

        if !allele_written
            write_allele_file(writers.allele_fh, all_vars, seq_id, location, ctx.sample_id_map)
            allele_written = true
        end
        write_transcript_product(writers.transcript_product_fh, all_vars, annotation, seq_id, location, ctx.sample_id_map)

        # Collect per-sample CANN entry keyed by original VCF alt allele (not IUPAC-derived base).
        # v.base may be an IUPAC ambiguity code for het calls; the VCF output write loop below
        # iterates record.alts, so we must use the original allele strings as keys here.
        for (alt, smap) in collect_cann_entries_for_annotation(variations, annotation, record, ctx.reference_strain, all_strains, strain_idx_map)
            for (strain, entries) in smap
                strain_map = get!(alt_strain_entries, alt, Dict{String, Vector{String}}())
                append!(get!(strain_map, strain, String[]), entries)
            end
        end
    end

    any_output || return false

    is_coding = any(a.is_coding == 1 for a in annotations) ? 1 : 0
    write_snp_feature(writers.snp_fh, first_all_vars, is_coding, seq_id, location,
                      ctx.reference_strain, ctx.all_strains)

    # HSSS binary files are SNP-only; skip positions where no record has a SNP-length allele.
    if any(has_snp_allele, records)
        unique_prods   = unique(all_annotation_products)
        hsss_prod_code = length(unique_prods) == 1 ?
            Int8(codepoint(first(only(unique_prods)))) : Int8(0)
        write_hsss_position!(writers.hsss, first_all_vars, ctx.reference_strain,
                             seq_id, location, ctx.all_strains, hsss_prod_code)  # ctx.all_strains is non-ref only; ref handled via ref_vars inside write_hsss_position!
    end

    (ref_keys, ref_cann_entries) = build_ref_cann_entries(annotations)
    modified_sample_data = fill_missing_coverage_gt(record, all_strains, chrom_coverage, ctx.ploidy)
    (alt_cann_entries, alt_strain_to_ca) = assign_cann_keys(alt_strain_entries, all_strains)

    # Write one VCF output entry per unique alt (for SnpEff downstream).
    n_orig_alts = length(record.alts)
    for (alt_i, alt) in enumerate(record.alts)
        if alt == "*"
            @warn "Unexpected * allele in VCF output at $(seq_id):$(location) — should have been removed by mergeVariantsByLocation.py"
            continue
        end
        haskey(alt_cann_entries, alt) || continue

        coding_alt_entries = filter(!=((".")), alt_cann_entries[alt])
        all_entries = [ref_cann_entries..., coding_alt_entries...]
        full_cann   = isempty(all_entries) ? "." : join(all_entries, ',')

        strain_to_alt_keys = get(alt_strain_to_ca, alt, Dict{String, Vector{String}}())
        ca_values = build_ca_values(record.format_keys, modified_sample_data, record.alts, alt,
                                    ref_keys, all_strains, strain_to_alt_keys)
        strain_to_dfs = Dict{String,String}(v.strain => string(v.downstream_of_frameshift) for v in variations)
        dfs_values = [get(strain_to_dfs, s, ".") for s in all_strains]
        write_vcf_entry(writers.vcf_fh, seq_id, location, record.ref, alt,
                        full_cann, record.info, record.format_keys,
                        modified_sample_data, ca_values, dfs_values,
                        n_orig_alts, alt_i)
    end

    true
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    args = parse_args(ARGS)
    global DEBUG = haskey(args, "debug")
    debug_log("Debug mode enabled")

    # Open VCF and parse header
    debug_log("Opening VCF: ", args["vcf_file"])
    (vcf_pf, all_strains, chrom_rank, info_headers) = open_vcf_peeked(args["vcf_file"])
    debug_log("VCF: ", length(all_strains), " strains")
    write_sample_dat([all_strains..., args["reference_strain"]])

    # Open VCF cache (may be absent/empty on first run)
    debug_log("Opening cache: ", args["cache_file"])
    cache_pf = open_cache_peeked(args["cache_file"])

    # Open coverage file
    debug_log("Opening coverage: ", args["coverage_file"])
    coverage_fh = open_coverage_file(args["coverage_file"])

    # Initialize processing context
    ctx = initialize_processing_context(args, all_strains)
    debug_log("Context: ", length(ctx.all_strains), " strains, ",
              length(ctx.cds_intervals), " CDS intervals")

    # Open output writers and write VCF header
    # Build full ordered strain list: reference first, then VCF non-ref strains in order
    all_strains_with_ref = [args["reference_strain"], all_strains...]
    writers = open_output_writers(args["output_vcf"], args["reference_strain"], all_strains_with_ref)
    write_vcf_cache_header(writers.vcf_fh, ctx.all_strains, info_headers)

    transcript_cache = TranscriptSequenceCache(Dict{String, Dict{String,String}}())

    chrom_coverage = Dict{String, Vector{Tuple{Int, Int, Float64}}}()
    current_chrom  = ""

    n_processed = 0

    while !vcf_pf.exhausted
        vcf_start = peek_sort_key(vcf_pf.line, chrom_rank)
        cache_key = cache_pf.exhausted ? (typemax(Int), typemax(Int)) :
                                          peek_sort_key(cache_pf.line, chrom_rank)

        # Drain cache entries that precede the current VCF record start
        # (positions that were variant in a prior run but are now absent)
        if !cache_pf.exhausted && cache_key < vcf_start
            advance!(cache_pf)
            continue
        end

        # Parse and advance VCF — collect all records sharing this chrom+pos
        records = VCFRecord[]
        while !vcf_pf.exhausted && peek_sort_key(vcf_pf.line, chrom_rank) == vcf_start
            push!(records, parse_vcf_record(vcf_pf.line, length(all_strains)))
            advance!(vcf_pf)
        end

        # Load coverage intervals when the chromosome changes
        if records[1].chrom != current_chrom
            current_chrom = records[1].chrom
            load_chrom_coverage!(coverage_fh, current_chrom, chrom_rank, chrom_coverage)
        end

        # Collect all cache entries at this (chrom, pos)
        cache_entries = Tuple{String,Int}[]
        while !cache_pf.exhausted
            ck = peek_sort_key(cache_pf.line, chrom_rank)
            ck != vcf_start && break
            parsed = parse_cache_tsv_record(cache_pf.line)
            if !isnothing(parsed)
                (_, _, tid, pos_in_cds) = parsed
                push!(cache_entries, (tid, pos_in_cds))
            end
            advance!(cache_pf)
        end

        if handle_variant_record!(records, cache_entries, ctx, writers, transcript_cache, all_strains, chrom_coverage)
            n_processed += 1
            if n_processed % 1000 == 0
                @info "Processed $n_processed variant positions"
            end
        end
    end

    debug_log("Processing complete. Total positions processed: ", n_processed)

    close_peeked(vcf_pf)
    close_peeked(cache_pf)
    close(coverage_fh.fh)
    close_output_writers(writers)
    close_processing_context(ctx)

    debug_log("Done!")
end

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
