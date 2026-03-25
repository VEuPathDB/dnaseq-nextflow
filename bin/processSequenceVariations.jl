#!/usr/bin/env julia

# processSequenceVariations.jl
# Reads a merged multi-sample FreeBayes GVCF and a coordinate-sorted VCF cache file,
# streams them concurrently in a sorted merge, annotates coding variants via SQLite
# transcript/indel databases, and writes four output files:
#   cache.vcf (CANN-annotated VCF cache), snpFeature.dat, allele.dat, product.dat

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

function find_cds(intervals::Vector{CDSInterval}, seq_id::String, location::Int)
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

    candidate == 0 && return nothing

    iv = intervals[candidate]
    if iv.seq_id == seq_id && location >= iv.start_pos && location <= iv.end_pos
        return iv
    end
    nothing
end

function compute_position_in_cds(tinfo::TranscriptInfo, location::Int, cds_number::Int)
    prior_len = 0
    for (idx, exon) in enumerate(tinfo.exons)
        if idx == cds_number
            return prior_len + (location - exon.start_pos) + 1
        end
        prior_len += exon.end_pos - exon.start_pos
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
# GVCF record structure
# ---------------------------------------------------------------------------

struct GVCFRecord
    chrom::String
    pos::Int
    ref::String
    alts::Vector{String}
    is_ref_block::Bool
    end_pos::Int                  # = pos for variants; = INFO/END for REF blocks
    info::String
    format_keys::Vector{String}
    sample_data::Vector{String}   # raw per-sample FORMAT strings
end

# ---------------------------------------------------------------------------
# VCF cache entry (result of parsing one cache data line)
# ---------------------------------------------------------------------------

struct CacheEntry
    cann_str::String
    ref_codon::String
    cds_number::Int
end

# ---------------------------------------------------------------------------
# Variation record (mutable, used during processing)
# ---------------------------------------------------------------------------

mutable struct Variation
    sequence_source_id::String
    location::Int
    strain::String
    reference::String
    base::String
    coverage::String
    percent::String
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
end

function Variation()
    Variation("", 0, "", "", "", "", "", "", "", "",
              0, 0, 0, 0, "", String[], "", "", 0, 0, 0)
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
    min_coverage::Int
end

struct OutputWriters
    vcf_cache_fh::IO
    snp_fh::IO
    allele_fh::IO
    product_fh::IO
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
    current_transcript_id::String
    current_transcript_seqs::Dict{String,String}
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
# GVCF I/O
# ---------------------------------------------------------------------------

"""
    parse_gvcf_header(io) -> (all_strains, chrom_rank)

Reads ## meta lines, builds chrom_rank from ##contig lines, extracts sample
names from #CHROM line. Leaves io positioned at first data line.
"""
function parse_gvcf_header(io::IO)
    chrom_rank = Dict{String,Int}()
    all_strains = String[]
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
            # Columns 10+ (1-based) are sample names
            all_strains = String[String(fields[i]) for i in 10:length(fields)]
            break
        end
    end

    debug_log("GVCF header: ", length(all_strains), " samples, ",
              length(chrom_rank), " contigs")
    (all_strains, chrom_rank, info_headers)
end

"""
    open_gvcf_peeked(path) -> (PeekedFile, all_strains, chrom_rank, info_headers)

Opens a bgzip-compressed GVCF via subprocess, parses its header, returns
a PeekedFile positioned at the first data line.
"""
function open_gvcf_peeked(path::String)
    io = open(`bgzip -d -c $path`)
    (all_strains, chrom_rank, info_headers) = parse_gvcf_header(io)
    pf = PeekedFile(io, "", false)
    advance!(pf)
    (pf, all_strains, chrom_rank, info_headers)
end

"""
    parse_gvcf_record(line, n_samples) -> GVCFRecord
"""
function parse_gvcf_record(line::String, n_samples::Int)::GVCFRecord
    fields = split(line, '\t')
    chrom = String(fields[1])
    pos   = parse(Int, fields[2])
    ref   = String(fields[4])
    alts  = String[String(a) for a in split(fields[5], ',')]
    info  = String(fields[8])
    fmt   = String(fields[9])

    is_ref_block = all(startswith(a, "<") for a in alts)

    end_pos = pos
    if is_ref_block
        m = match(r"END=(\d+)", info)
        !isnothing(m) && (end_pos = parse(Int, m.captures[1]))
    end

    format_keys = String[String(k) for k in split(fmt, ':')]
    sample_data = String[String(fields[9+i]) for i in 1:n_samples if 9+i <= length(fields)]

    GVCFRecord(chrom, pos, ref, alts, is_ref_block, end_pos, info, format_keys, sample_data)
end

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
# VCF cache I/O
# ---------------------------------------------------------------------------

"""
    open_cache_peeked(path) -> PeekedFile

Opens the VCF cache file for reading. Returns a pre-exhausted PeekedFile
if the file is absent or empty. Skips VCF header lines (#-prefixed).
"""
function open_cache_peeked(path::String)::PeekedFile
    if !isfile(path) || filesize(path) == 0
        return PeekedFile(IOBuffer(""), "", true)
    end

    fh = open(path, "r")
    # Skip header lines
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
    # File contained only headers
    close(fh)
    PeekedFile(IOBuffer(""), "", true)
end

"""
    parse_cache_vcf_record(line) -> (chrom, pos, ref, alt, cann_str, ref_codon, cds_number) or nothing
"""
function parse_cache_vcf_record(line::String)
    startswith(line, '#') && return nothing
    fields = split(line, '\t')
    length(fields) < 8 && return nothing

    chrom = String(fields[1])
    pos   = parse(Int, fields[2])
    ref   = String(fields[4])
    alt   = String(fields[5])
    info  = String(fields[8])

    cann_str   = "."
    ref_codon  = ""
    cds_number = 0

    for part in split(info, ';')
        if startswith(part, "CANN=")
            cann_str = String(part[6:end])
        elseif startswith(part, "RC=")
            ref_codon = String(part[4:end])
        elseif startswith(part, "CDSNR=")
            cds_number = parse(Int, part[7:end])
        end
    end

    (chrom, pos, ref, alt, cann_str, ref_codon, cds_number)
end

function write_vcf_cache_header(fh::IO, all_strains::Vector{String}, info_headers::Vector{String})
    write(fh, "##fileformat=VCFv4.2\n")
    for h in info_headers
        write(fh, h, "\n")
    end
    write(fh, "##INFO=<ID=CANN,Number=.,Type=String,Description=\"Coding annotation: key:codon:alt_aa:effect:transcript_id:pos_in_cds:pos_in_codon\">\n")
    write(fh, "##INFO=<ID=RC,Number=1,Type=String,Description=\"Reference codon\">\n")
    write(fh, "##INFO=<ID=CDSNR,Number=1,Type=Integer,Description=\"CDS exon number\">\n")
    chrom_line = join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", all_strains...], '\t')
    write(fh, chrom_line, "\n")
end

function write_vcf_cache_entry(fh::IO, chrom::String, pos::Int, ref::String, alt::String,
                                cann_str::String, ref_codon::String, cds_number::Int,
                                gvcf_info::String, format_keys::Vector{String}, sample_data::Vector{String})
    info   = "$(gvcf_info);CANN=$(cann_str);RC=$(ref_codon);CDSNR=$(cds_number)"
    format = join(format_keys, ':')
    write(fh, join([chrom, string(pos), ".", ref, alt, ".", ".", info, format, sample_data...], '\t'), "\n")
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
    gt_allele_idx(gt) -> Int

Returns the primary alt allele index from GT (0 = ref). Used for percent computation.
"""
function gt_allele_idx(gt::String)::Int
    sep_idx = findfirst(c -> c == '/' || c == '|', gt)
    isnothing(sep_idx) && return parse(Int, gt)
    parse(Int, gt[1:sep_idx-1])
end

"""
    compute_percent(fmt, allele_idx) -> String

Computes AO/(RO+AO)*100 for alt allele. Returns "0.0" for ref or missing.
"""
function compute_percent(fmt::Dict{String,String}, allele_idx::Int)::String
    allele_idx == 0 && return "0.0"

    ao_str = get(fmt, "AO", "")
    ro_str = get(fmt, "RO", "0")
    isempty(ao_str) && return "0.0"

    ao_values = split(ao_str, ',')
    allele_idx > length(ao_values) && return "0.0"

    ao = parse(Float64, ao_values[allele_idx])
    ro = parse(Float64, isempty(ro_str) ? "0" : ro_str)
    total = ro + ao
    total <= 0 && return "0.0"

    @sprintf("%.2f", ao / total * 100)
end

"""
    build_variations_from_record(record, all_strains, undone_strains, min_coverage)
        -> Vector{Variation}

Builds per-strain Variation records from a GVCF variant record.
Skips undone strains, missing GTs, and samples below min_coverage.
"""
function build_variations_from_record(
    record::GVCFRecord,
    all_strains::Vector{String},
    undone_strains::Set{String},
    min_coverage::Int
)::Vector{Variation}
    variations = Variation[]

    for (i, strain) in enumerate(all_strains)
        strain in undone_strains && continue
        i > length(record.sample_data) && continue

        fmt = parse_format_field(record.format_keys, record.sample_data[i])

        gt = get(fmt, "GT", "")
        (isempty(gt) || gt == "." || gt == "./." || gt == ".|.") && continue

        dp_str = get(fmt, "DP", "0")
        dp = isempty(dp_str) || dp_str == "." ? 0 : parse(Int, dp_str)
        dp < min_coverage && continue

        base = gt_to_base(gt, record.ref, record.alts)
        isempty(base) && continue

        aidx = gt_allele_idx(gt)
        pct  = compute_percent(fmt, aidx)
        gq   = get(fmt, "GQ", "0")

        v = Variation()
        v.sequence_source_id = record.chrom
        v.location           = record.pos
        v.strain             = strain
        v.reference          = record.ref
        v.base               = base
        v.coverage           = string(dp)
        v.percent            = pct
        v.quality            = gq
        v.pvalue             = "."
        v.snp_source_id      = "NGS_SNP.$(record.chrom).$(record.pos)"
        v.matches_reference  = (base == record.ref) ? 1 : 0

        push!(variations, v)
    end

    variations
end

# ---------------------------------------------------------------------------
# Sorted-merge helpers: sort keys over VCF/GVCF lines
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

"""
    peek_end_key(line, chrom_rank) -> (Int, Int)

Returns (chrom_rank, END) for GVCF REF blocks; (chrom_rank, POS) for variants.
Used to determine the span of a REF block without full parsing.
"""
function peek_end_key(line::String, chrom_rank::Dict{String,Int})::Tuple{Int,Int}
    fields = split(line, '\t'; limit=9)
    length(fields) < 8 && return peek_sort_key(line, chrom_rank)

    chrom   = String(fields[1])
    pos     = parse(Int, fields[2])
    alt_str = String(fields[5])
    info    = String(fields[8])
    rank    = get(chrom_rank, chrom, typemax(Int))

    if startswith(alt_str, "<")
        m = match(r"END=(\d+)", info)
        !isnothing(m) && return (rank, parse(Int, m.captures[1]))
    end

    (rank, pos)
end

# ---------------------------------------------------------------------------
# Resource management
# ---------------------------------------------------------------------------

"""
    initialize_processing_context(args, all_strains, min_coverage) -> ProcessingContext
"""
function initialize_processing_context(args, all_strains::Vector{String}, min_coverage::Int)
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

    ProcessingContext(
        args["reference_strain"],
        undone_strains,
        cds_intervals,
        transcript_info,
        transcript_db,
        indel_db,
        fs_info,
        all_strains,
        min_coverage
    )
end

"""
    open_output_writers(output_vcf) -> OutputWriters
"""
function open_output_writers(output_vcf)
    vcf_cache_fh = open(output_vcf, "w")
    snp_fh       = open("snpFeature.dat", "w")
    allele_fh    = open("allele.dat", "w")
    product_fh   = open("product.dat", "w")

    OutputWriters(vcf_cache_fh, snp_fh, allele_fh, product_fh)
end

function close_processing_context(ctx::ProcessingContext)
    close(ctx.transcript_db)
    close(ctx.indel_db)
end

function close_output_writers(writers::OutputWriters)
    close(writers.vcf_cache_fh)
    close(writers.snp_fh)
    close(writers.allele_fh)
    close(writers.product_fh)
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
    annotate_position(seq_id, location, ctx, transcript_cache) -> PositionAnnotation

CDS lookup and codon/position computation. Updates transcript_cache if transcript changes.
"""
function annotate_position(
    seq_id::String,
    location::Int,
    ctx::ProcessingContext,
    transcript_cache::TranscriptSequenceCache
)
    cds_hit = find_cds(ctx.cds_intervals, seq_id, location)

    is_coding     = 0
    transcript_id = ""
    cds_number    = 0
    pos_in_cds    = 0
    pos_in_codon_val = 0
    ref_codon     = ""
    ref_product   = ""

    if !isnothing(cds_hit)
        is_coding     = 1
        transcript_id = cds_hit.transcript_id
        cds_number    = cds_hit.cds_number

        tinfo = ctx.transcript_info[transcript_id]
        pos_in_cds       = compute_position_in_cds(tinfo, location, cds_number)
        pos_in_codon_val = position_in_codon(pos_in_cds)

        if transcript_id != transcript_cache.current_transcript_id
            debug_log("    Loading sequences for transcript: ", transcript_id)
            transcript_cache.current_transcript_id   = transcript_id
            transcript_cache.current_transcript_seqs = load_transcript_sequences(ctx.transcript_db, transcript_id)
        end

        ref_seq = get(transcript_cache.current_transcript_seqs, ctx.reference_strain, "")
        if !isempty(ref_seq)
            ref_codon   = extract_codon(ref_seq, pos_in_cds)
            ref_product = translate_codon(ref_codon)
        end
    end

    PositionAnnotation(is_coding, transcript_id, cds_number, pos_in_cds,
                       pos_in_codon_val, ref_codon, ref_product)
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

        strain_seq = get(transcript_cache.current_transcript_seqs, v.strain, "")
        if !isempty(strain_seq)
            shift        = get_indel_shift(ctx.indel_db, annotation.transcript_id, v.strain, annotation.pos_in_cds)
            adjusted_pos = annotation.pos_in_cds + shift
            strain_codon = extract_codon(strain_seq, adjusted_pos)
            v.codon   = strain_codon
            expanded  = expand_codon(strain_codon)
            v.product = [translate_codon(c) for c in expanded]
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
        1
    )
end

# ---------------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------------

function write_snp_feature(
    snp_fh::IO,
    variations::Vector{Variation},
    annotation::PositionAnnotation,
    seq_id::String,
    location::Int,
    reference_strain::String
)
    ref_allele = variations[1].reference
    isempty(ref_allele) && (ref_allele = variations[1].base)

    allele_counts   = Dict{String,Int}()
    product_counts  = Dict{String,Int}()
    strain_set      = Set{String}()
    has_stop_codon  = 0
    total_allele_count = length(variations)

    for v in variations
        allele_counts[v.base] = get(allele_counts, v.base, 0) + 1
        push!(strain_set, v.strain)
        for p in v.product
            product_counts[p] = get(product_counts, p, 0) + 1
            p == "*" && (has_stop_codon = 1)
        end
    end

    distinct_strain_count  = length(strain_set)
    distinct_allele_count  = length(allele_counts)
    has_nonsynonymous_allele = length(product_counts) > 1 ? 1 : 0

    sorted_alleles  = sort(collect(keys(allele_counts));  by = a -> (-allele_counts[a], a))
    sorted_products = sort(collect(keys(product_counts)); by = p -> (-product_counts[p], p))

    major_allele       = sorted_alleles[1]
    minor_allele       = length(sorted_alleles) > 1 ? sorted_alleles[2] : ""
    major_allele_count = allele_counts[major_allele]
    minor_allele_count = length(sorted_alleles) > 1 ? allele_counts[minor_allele] : ""

    major_product = length(sorted_products) > 0 ? sorted_products[1] : ""
    minor_product = length(sorted_products) > 1 ? sorted_products[2] : ""

    write(snp_fh, join([
        string(location),
        annotation.transcript_id,
        seq_id,
        reference_strain,
        ref_allele,
        string(has_nonsynonymous_allele),
        major_allele,
        minor_allele,
        string(major_allele_count),
        string(minor_allele_count),
        major_product,
        minor_product,
        string(distinct_strain_count),
        string(distinct_allele_count),
        string(annotation.is_coding),
        string(total_allele_count),
        string(has_stop_codon),
        annotation.ref_codon
    ], "\t"), "\n")
end

function write_allele_and_product_files(
    allele_fh::IO,
    product_fh::IO,
    variations::Vector{Variation},
    annotation::PositionAnnotation
)
    annotation.is_coding != 1 && return

    allele_counts = Dict{String,Int}()
    for v in variations
        allele_counts[v.base] = get(allele_counts, v.base, 0) + 1
    end

    for allele in keys(allele_counts)
        distinct_strains = Set{String}()
        allele_count = 0
        sum_coverage = 0.0
        sum_percent  = 0.0

        for v in variations
            v.base != allele && continue
            allele_count += 1
            push!(distinct_strains, v.strain)
            cov = isempty(v.coverage) ? 0.0 : parse(Float64, v.coverage)
            pct = isempty(v.percent)  ? 0.0 : parse(Float64, v.percent)
            sum_coverage += cov
            sum_percent  += pct
        end

        avg_cov = allele_count > 0 ? sum_coverage / allele_count : 0.0
        avg_pct = allele_count > 0 ? sum_percent  / allele_count : 0.0

        write(allele_fh, join([
            allele,
            string(length(distinct_strains)),
            string(allele_count),
            @sprintf("%.2f", avg_cov),
            @sprintf("%.2f", avg_pct)
        ], "\t"), "\n")
    end

    all_product_counts = Dict{String,Int}()
    for v in variations
        for p in v.product
            all_product_counts[p] = get(all_product_counts, p, 0) + 1
        end
    end

    for v in variations
        isempty(v.codon) && continue
        expanded_codons = expand_codon(v.codon)
        for ec in expanded_codons
            product = translate_codon(ec)
            count   = get(all_product_counts, product, 0)
            write(product_fh, join([
                ec,
                string(annotation.pos_in_codon_val),
                annotation.transcript_id,
                string(count),
                product,
                string(annotation.pos_in_cds),
                string(annotation.pos_in_codon_val),
                string(v.downstream_of_frameshift)
            ], "\t"), "\n")
        end
    end
end

# ---------------------------------------------------------------------------
# CANN annotation
# ---------------------------------------------------------------------------

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
    is_indel   = (ref_len != alt_len)
    is_complex = is_indel && !isempty(ref_allele) && !isempty(alt_allele) && ref_allele[1] != alt_allele[1]
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
        return "k0:.:.:$(structural):$(tid):$(pos_in_cds):$(pic)"
    end

    # SNP or complex variant: compute amino acid effect
    codon       = isempty(v.codon) ? "." : v.codon
    unique_prods = unique(v.product)
    product_str  = isempty(unique_prods) ? "." : join(unique_prods, "/")

    has_stop = any(p == "*" for p in unique_prods)
    aa_effect = if has_stop
        "nonsense"
    elseif product_str != annotation.ref_product
        "missense"
    else
        "synonymous"
    end

    if !is_indel
        return "k0:$(codon):$(product_str):$(aa_effect):$(tid):$(pos_in_cds):$(pic)"
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
        return "k0:$(codon):$(product_str):$(aa_effect);$(structural):$(tid):$(pos_in_cds):$(pic)"
    end
end

"""
    decode_cann_to_annotation(ce, ctx, transcript_cache) -> PositionAnnotation

Reconstructs a PositionAnnotation from a cached CacheEntry.
Also loads transcript sequences into transcript_cache so annotate_variations!
can access them.
"""
function decode_cann_to_annotation(
    ce::CacheEntry,
    ctx::ProcessingContext,
    transcript_cache::TranscriptSequenceCache
)::PositionAnnotation
    ce.cann_str == "." && return PositionAnnotation(0, "", 0, 0, 0, "", "")

    entry = split(ce.cann_str, ',')[1]
    parts = split(entry, ':')
    # parts: key, codon, alt_aa, effect, transcript_id, pos_in_cds, pos_in_codon
    length(parts) < 7 && return PositionAnnotation(0, "", 0, 0, 0, "", "")

    transcript_id    = String(parts[5])
    pos_in_cds       = parts[6] == "." ? 0 : parse(Int, parts[6])
    pos_in_codon_val = parts[7] == "." ? 0 : parse(Int, parts[7])
    ref_codon        = ce.ref_codon
    ref_product      = isempty(ref_codon) ? "" : translate_codon(ref_codon)

    if !isempty(transcript_id) && transcript_id != transcript_cache.current_transcript_id
        debug_log("    Loading sequences for transcript (cache hit): ", transcript_id)
        transcript_cache.current_transcript_id   = transcript_id
        transcript_cache.current_transcript_seqs = load_transcript_sequences(ctx.transcript_db, transcript_id)
    end

    PositionAnnotation(1, transcript_id, ce.cds_number, pos_in_cds,
                       pos_in_codon_val, ref_codon, ref_product)
end

# ---------------------------------------------------------------------------
# Core record handler
# ---------------------------------------------------------------------------

"""
    handle_variant_record!(record, cache_entries, ctx, writers, transcript_cache, all_strains) -> Bool

Processes one variant GVCF record end-to-end. Returns true if output was written.
cache_entries: Dict keyed by (ref, alt) for positions with a cache hit at this coordinate.
"""
function handle_variant_record!(
    record::GVCFRecord,
    cache_entries::Dict{Tuple{String,String},CacheEntry},
    ctx::ProcessingContext,
    writers::OutputWriters,
    transcript_cache::TranscriptSequenceCache,
    all_strains::Vector{String}
)::Bool
    seq_id   = record.chrom
    location = record.pos
    debug_log("Processing position: ", seq_id, ":", location)

    variations = build_variations_from_record(record, all_strains, ctx.undone_strains, ctx.min_coverage)
    isempty(variations) && return false

    # Determine position annotation: cache hit or fresh CDS lookup
    annotation = if !isempty(cache_entries)
        ce = first(values(cache_entries))
        decode_cann_to_annotation(ce, ctx, transcript_cache)
    else
        annotate_position(seq_id, location, ctx, transcript_cache)
    end

    annotate_variations!(variations, annotation, ctx, transcript_cache)

    # Set downstream_of_frameshift on each variation
    for v in variations
        if !isempty(v.transcript)
            v.downstream_of_frameshift = check_downstream_of_frameshift(
                ctx.fs_info, v.strain, v.transcript, v.location)
        end
    end

    ref_variation = build_reference_variation(variations, annotation, seq_id, location, ctx.reference_strain)
    push!(variations, ref_variation)

    !has_variation(variations) && return false

    # Write VCF cache entries: one per unique alt allele seen in this record
    seen_alts = Set{String}()
    for v in variations
        v.strain == ctx.reference_strain && continue
        v.base == record.ref && continue  # ref call
        alt = v.base
        alt in seen_alts && continue
        push!(seen_alts, alt)

        cann_str = build_cann_string(record.ref, alt, v, annotation)
        write_vcf_cache_entry(writers.vcf_cache_fh, seq_id, location, record.ref, alt,
                              cann_str, annotation.ref_codon, annotation.cds_number,
                              record.info, record.format_keys, record.sample_data)
    end

    write_snp_feature(writers.snp_fh, variations, annotation, seq_id, location, ctx.reference_strain)
    write_allele_and_product_files(writers.allele_fh, writers.product_fh, variations, annotation)

    true
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    args = parse_args(ARGS)
    global DEBUG = haskey(args, "debug")
    debug_log("Debug mode enabled")

    min_coverage = haskey(args, "min_coverage") ? parse(Int, args["min_coverage"]) : 1

    # Open GVCF and parse header
    debug_log("Opening GVCF: ", args["vcf_file"])
    (gvcf_pf, all_strains, chrom_rank, info_headers) = open_gvcf_peeked(args["vcf_file"])
    debug_log("GVCF: ", length(all_strains), " strains")

    # Open VCF cache (may be absent/empty on first run)
    debug_log("Opening cache: ", args["cache_file"])
    cache_pf = open_cache_peeked(args["cache_file"])

    # Initialize processing context
    ctx = initialize_processing_context(args, all_strains, min_coverage)
    debug_log("Context: ", length(ctx.all_strains), " strains, ",
              length(ctx.cds_intervals), " CDS intervals")

    # Open output writers and write VCF cache header
    writers = open_output_writers(args["output_vcf"])
    write_vcf_cache_header(writers.vcf_cache_fh, ctx.all_strains, info_headers)

    transcript_cache = TranscriptSequenceCache("", Dict{String,String}())

    # Span-aware sorted-merge loop
    n_processed = 0

    while !gvcf_pf.exhausted
        gvcf_start = peek_sort_key(gvcf_pf.line, chrom_rank)
        gvcf_end   = peek_end_key(gvcf_pf.line, chrom_rank)
        cache_key  = cache_pf.exhausted ? (typemax(Int), typemax(Int)) :
                                           peek_sort_key(cache_pf.line, chrom_rank)

        # Drain cache entries that precede the current GVCF record start
        # (positions that were variant in a prior run but are now absent)
        if !cache_pf.exhausted && cache_key < gvcf_start
            advance!(cache_pf)
            continue
        end

        # Parse and advance GVCF
        record = parse_gvcf_record(gvcf_pf.line, length(all_strains))
        advance!(gvcf_pf)

        if record.is_ref_block
            # Drain all cache entries within this REF block span [pos, end_pos]
            # These positions were variant before but are now reference-covered
            while !cache_pf.exhausted
                ck = peek_sort_key(cache_pf.line, chrom_rank)
                ck > gvcf_end && break
                advance!(cache_pf)
            end
            continue
        end

        # Variant record: collect all cache entries at this (chrom, pos)
        cache_entries = Dict{Tuple{String,String},CacheEntry}()
        while !cache_pf.exhausted
            ck = peek_sort_key(cache_pf.line, chrom_rank)
            ck != gvcf_start && break
            parsed = parse_cache_vcf_record(cache_pf.line)
            if !isnothing(parsed)
                (_, _, ref, alt, cann_str, ref_codon, cds_number) = parsed
                cache_entries[(ref, alt)] = CacheEntry(cann_str, ref_codon, cds_number)
            end
            advance!(cache_pf)
        end

        if handle_variant_record!(record, cache_entries, ctx, writers, transcript_cache, all_strains)
            n_processed += 1
            if n_processed % 1000 == 0
                @info "Processed $n_processed variant positions"
            end
        end
    end

    debug_log("Processing complete. Total positions processed: ", n_processed)

    close_peeked(gvcf_pf)
    close_peeked(cache_pf)
    close_output_writers(writers)
    close_processing_context(ctx)

    debug_log("Done!")
end

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
main()
