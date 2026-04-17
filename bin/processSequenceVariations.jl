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

    fh = endswith(path, ".gz") ? open(`bgzip -d -c $path`) : open(path, "r")
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
    parse_cache_vcf_record(line) -> (chrom, pos, ref, alt, cann_str) or nothing
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

    cann_str = "."
    for part in split(info, ';')
        if startswith(part, "CANN=")
            cann_str = String(part[6:end])
            break
        end
    end

    (chrom, pos, ref, alt, cann_str)
end

function write_vcf_cache_header(fh::IO, all_strains::Vector{String}, info_headers::Vector{String})
    write(fh, "##fileformat=VCFv4.2\n")
    for h in info_headers
        write(fh, h, "\n")
    end
    write(fh, "##INFO=<ID=CANN,Number=.,Type=String,Description=\"Coding annotation entries, comma-separated. r-prefixed keys (r0,r1,...) = reference allele per transcript; k-prefixed keys (k0,k1,...) = alt allele per transcript. Format per entry: key:codon:aa:effect:transcript_id:pos_in_cds:pos_in_codon. Compound effects use '&' separator (e.g. missense&frameshift).\">\n")
    write(fh, "##FORMAT=<ID=CA,Number=1,Type=String,Description=\"CANN key(s) per GT allele. Alleles separated by '/' (unphased) or '|' (phased). Multiple transcript keys for one allele separated by ';'. 'r'=ref allele no CDS annotation, '.'=missing/no-call\">\n")
    write(fh, "##FORMAT=<ID=DFS,Number=1,Type=Integer,Description=\"Downstream of frameshift: 1 if this sample carries an upstream indel that disrupts the reading frame at this position, 0 otherwise.\">\n")
    chrom_line = join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", all_strains...], '\t')
    write(fh, chrom_line, "\n")
end

function write_vcf_cache_entry(fh::IO, chrom::String, pos::Int, ref::String, alt::String,
                                cann_str::String, gvcf_info::String,
                                format_keys::Vector{String}, sample_data::Vector{String},
                                ca_values::Vector{String}, dfs_values::Vector{String})
    info    = "$(gvcf_info);CANN=$(cann_str)"
    format  = join([format_keys..., "CA", "DFS"], ':')
    samples = [i <= length(ca_values) ? "$(sample_data[i]):$(ca_values[i]):$(get(dfs_values, i, "."))" : "$(sample_data[i]):.:."
               for i in 1:length(sample_data)]
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
    build_variations_from_record(record, all_strains, undone_strains, prev_coverage_span)
        -> Vector{Variation}

Builds per-strain Variation records from a GVCF variant record.
Skips undone strains and missing GTs.
"""
function build_variations_from_record(
    record::GVCFRecord,
    all_strains::Vector{String},
    undone_strains::Set{String},
    prev_coverage_span::Dict{String, Tuple{String, Int, Int}}
)::Vector{Variation}
    variations = Variation[]

    for (i, strain) in enumerate(all_strains)
        strain in undone_strains && continue
        i > length(record.sample_data) && continue

        fmt = parse_format_field(record.format_keys, record.sample_data[i])

        gt = get(fmt, "GT", "")
        if isempty(gt) || gt == "." || gt == "./." || gt == ".|."
            # Sample has no call at this position; check if a prior REF block covers it
            span = get(prev_coverage_span, strain, nothing)
            if !isnothing(span)
                span_chrom, span_end, span_dp = span
                if span_chrom == record.chrom && record.pos <= span_end
                    v = Variation()
                    v.sequence_source_id = record.chrom
                    v.location           = record.pos
                    v.strain             = strain
                    v.reference          = record.ref
                    v.base               = record.ref
                    v.coverage           = string(span_dp)
                    v.percent            = "100"
                    v.quality            = "."
                    v.pvalue             = "."
                    v.snp_source_id      = "NGS_SNP.$(record.chrom).$(record.pos)"
                    v.matches_reference  = 1
                    push!(variations, v)
                end
            end
            continue
        end

        dp_str = get(fmt, "DP", "0")
        dp     = isempty(dp_str) || dp_str == "." ? 0 : parse(Int, dp_str)

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

    ProcessingContext(
        args["reference_strain"],
        undone_strains,
        cds_intervals,
        transcript_info,
        transcript_db,
        indel_db,
        fs_info,
        all_strains
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
                strain_codon = extract_codon(strain_seq, annotation.pos_in_cds)
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

# ---------------------------------------------------------------------------
# CANN annotation
# ---------------------------------------------------------------------------

"""
    build_ref_cann_entry(key, annotation) -> String

Builds the r-keyed CANN entry for the reference allele at a coding position.
Format mirrors alt entries: key:codon:aa:effect:transcript_id:pos_in_cds:pos_in_codon
"""
function build_ref_cann_entry(key::String, annotation::PositionAnnotation)::String
    annotation.is_coding != 1 && return "."
    codon = isempty(annotation.ref_codon)   ? "." : annotation.ref_codon
    aa    = isempty(annotation.ref_product) ? "." : annotation.ref_product
    "$(key):$(codon):$(aa):reference:$(annotation.transcript_id):$(annotation.pos_in_cds):$(annotation.pos_in_codon_val)"
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
        return "k0:.:.:$(structural):$(tid):$(pos_in_cds):$(pic)"
    end

    # SNP or complex variant: compute amino acid effect
    codon       = isempty(v.codon) ? "." : v.codon
    unique_prods = unique(v.product)
    product_str  = isempty(unique_prods) ? "." : join(unique_prods, "/")

    # Codon/product suppressed because strain is downstream of a frameshift
    if codon == "." && isempty(unique_prods)
        return "k0:.:.:downstream_frameshift:$(tid):$(pos_in_cds):$(pic)"
    end

    # Codon contains ambiguous base(s) — skip product and effect
    if occursin(r"[Nn]", codon)
        return "k0:$(codon):.:.:$(tid):$(pos_in_cds):$(pic)"
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
        return "k0:$(codon):$(product_str):$(aa_effect)&$(structural):$(tid):$(pos_in_cds):$(pic)"
    end
end

"""
    decode_all_cann_annotations(ce, ctx, transcript_cache) -> Vector{PositionAnnotation}

Reconstructs one PositionAnnotation per unique transcript encoded in ce.cann_str.
The comma-separated CANN entries may cover multiple overlapping transcripts.
Returns a single non-coding annotation when cann_str is ".".
"""
function decode_all_cann_annotations(
    ce::CacheEntry,
    ctx::ProcessingContext,
    transcript_cache::TranscriptSequenceCache
)::Vector{PositionAnnotation}
    ce.cann_str == "." && return [PositionAnnotation(0, "", 0, 0, 0, "", "")]

    # Reconstruct PositionAnnotations from r-keyed entries only.
    # r-keyed entries encode the reference allele annotation per transcript and carry
    # all the data needed (ref_codon, ref_product, pos_in_cds, pos_in_codon).
    # k-keyed entries (alt annotations) are ignored here — they will be recomputed
    # from SQLite sequences by annotate_variations!.
    seen_transcripts = Set{String}()
    annotations      = PositionAnnotation[]

    for entry in split(ce.cann_str, ',')
        parts = split(entry, ':')
        length(parts) < 7 && continue
        startswith(String(parts[1]), 'r') || continue   # only r-keyed entries

        transcript_id = String(parts[5])
        isempty(transcript_id) && continue
        transcript_id in seen_transcripts && continue
        push!(seen_transcripts, transcript_id)

        ref_codon        = parts[2] == "." ? "" : String(parts[2])
        ref_product      = parts[3] == "." ? "" : String(parts[3])
        pos_in_cds       = parts[6] == "." ? 0  : parse(Int, parts[6])
        pos_in_codon_val = parts[7] == "." ? 0  : parse(Int, parts[7])

        debug_log("    Loading sequences for transcript (cache hit): ", transcript_id)
        load_transcript_if_needed!(transcript_cache, ctx.transcript_db, transcript_id)

        push!(annotations, PositionAnnotation(1, transcript_id, 0, pos_in_cds,
                                              pos_in_codon_val, ref_codon, ref_product))
    end

    isempty(annotations) ? [PositionAnnotation(0, "", 0, 0, 0, "", "")] : annotations
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
    all_strains::Vector{String},
    prev_coverage_span::Dict{String, Tuple{String, Int, Int}}
)::Bool
    seq_id   = record.chrom
    location = record.pos
    debug_log("Processing position: ", seq_id, ":", location)

    variations = build_variations_from_record(record, all_strains, ctx.undone_strains, prev_coverage_span)
    isempty(variations) && return false

    # Determine annotations: one per overlapping transcript (cache or fresh GTF lookup)
    annotations = if !isempty(cache_entries)
        ce = first(values(cache_entries))
        decode_all_cann_annotations(ce, ctx, transcript_cache)
    else
        annotate_position_all(seq_id, location, ctx, transcript_cache)
    end

    # Accumulate per-sample CANN entries across all annotations (transcripts).
    # alt_strain_entries: alt -> strain -> [entry per transcript, in annotation order]
    alt_strain_entries = Dict{String, Dict{String, Vector{String}}}()
    first_annotation   = annotations[1]
    any_output         = false

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

        write_snp_feature(writers.snp_fh, all_vars, annotation, seq_id, location, ctx.reference_strain)
        write_allele_and_product_files(writers.allele_fh, writers.product_fh, all_vars, annotation)

        # Collect per-sample CANN entry for each alt under this annotation
        for v in variations
            v.strain == ctx.reference_strain && continue
            v.base == record.ref && continue
            entry = build_cann_string(record.ref, v.base, v, annotation)
            strain_map = get!(alt_strain_entries, v.base, Dict{String, Vector{String}}())
            push!(get!(strain_map, v.strain, String[]), entry)
        end
    end

    any_output || return false

    # Build r-keyed CANN entries for the reference allele (one per coding transcript).
    # These are identical across all alt cache lines at this position.
    ref_keys         = String[]
    ref_cann_entries = String[]
    for (i, annotation) in enumerate(annotations)
        annotation.is_coding == 0 && continue
        key   = "r$(length(ref_keys))"
        entry = build_ref_cann_entry(key, annotation)
        entry == "." && continue
        push!(ref_keys, key)
        push!(ref_cann_entries, entry)
    end

    # Build modified sample data: fill in GT=0 and DP for samples covered by a prior
    # REF block span but absent from this variant record (missing GT).
    gt_idx = findfirst(==("GT"), record.format_keys)
    dp_idx = findfirst(==("DP"), record.format_keys)
    modified_sample_data = copy(record.sample_data)
    for (i, strain) in enumerate(all_strains)
        i > length(modified_sample_data) && continue
        fmt = parse_format_field(record.format_keys, modified_sample_data[i])
        gt = get(fmt, "GT", "")
        (isempty(gt) || gt == "." || gt == "./." || gt == ".|.") || continue
        span = get(prev_coverage_span, strain, nothing)
        isnothing(span) && continue
        span_chrom, span_end, span_dp = span
        (span_chrom == record.chrom && record.pos <= span_end) || continue
        fields = fill(".", length(record.format_keys))
        !isnothing(gt_idx) && (fields[gt_idx] = "0")
        !isnothing(dp_idx) && (fields[dp_idx] = string(span_dp))
        modified_sample_data[i] = join(fields, ":")
    end

    # Build per-alt CANN key assignments and per-strain CA mappings.
    # Unique annotation entries across all strains are deduplicated; each unique
    # entry gets a k-key. Strains with the same annotation share a key.
    # k-keys are local to each alt (reset per alt).
    alt_cann_entries   = Dict{String, Vector{String}}()          # alt -> ordered unique keyed entries
    alt_strain_to_ca   = Dict{String, Dict{String, Vector{String}}}()  # alt -> strain -> [keys per transcript]

    for (alt, strain_map) in alt_strain_entries
        entry_to_key   = Dict{String, String}()
        ordered_keyed  = String[]

        for strain in all_strains  # canonical strain order for deterministic key assignment
            strain_entries = get(strain_map, strain, nothing)
            isnothing(strain_entries) && continue
            for entry in strain_entries
                entry == "." && continue
                if !haskey(entry_to_key, entry)
                    key = "k$(length(entry_to_key))"
                    entry_to_key[entry] = key
                    push!(ordered_keyed, replace(entry, r"^k0:" => "$(key):", count=1))
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

    # Write one VCF cache entry per unique alt.
    for alt in record.alts
        haskey(alt_cann_entries, alt) || continue

        coding_alt_entries = filter(!=((".")), alt_cann_entries[alt])
        all_entries = [ref_cann_entries..., coding_alt_entries...]
        full_cann   = isempty(all_entries) ? "." : join(all_entries, ',')

        strain_to_alt_keys = get(alt_strain_to_ca, alt, Dict{String, Vector{String}}())
        ca_values = build_ca_values(record.format_keys, modified_sample_data, record.alts, alt,
                                    ref_keys, all_strains, strain_to_alt_keys)
        strain_to_dfs = Dict{String,String}(v.strain => string(v.downstream_of_frameshift) for v in variations)
        dfs_values = [get(strain_to_dfs, s, ".") for s in all_strains]
        write_vcf_cache_entry(writers.vcf_cache_fh, seq_id, location, record.ref, alt,
                              full_cann, record.info, record.format_keys,
                              modified_sample_data, ca_values, dfs_values)
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

    # Open GVCF and parse header
    debug_log("Opening GVCF: ", args["vcf_file"])
    (gvcf_pf, all_strains, chrom_rank, info_headers) = open_gvcf_peeked(args["vcf_file"])
    debug_log("GVCF: ", length(all_strains), " strains")

    # Open VCF cache (may be absent/empty on first run)
    debug_log("Opening cache: ", args["cache_file"])
    cache_pf = open_cache_peeked(args["cache_file"])

    # Initialize processing context
    ctx = initialize_processing_context(args, all_strains)
    debug_log("Context: ", length(ctx.all_strains), " strains, ",
              length(ctx.cds_intervals), " CDS intervals")

    # Open output writers and write VCF cache header
    writers = open_output_writers(args["output_vcf"])
    write_vcf_cache_header(writers.vcf_cache_fh, ctx.all_strains, info_headers)

    transcript_cache = TranscriptSequenceCache(Dict{String, Dict{String,String}}())

    # Per-sample last-seen REF block span: strain -> (chrom, end_pos, dp)
    # Used to synthesize reference calls at variant positions for covered samples
    # that bcftools merge left as missing GT.
    prev_coverage_span = Dict{String, Tuple{String, Int, Int}}()

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
            # Update per-sample coverage spans from this REF block
            for (i, strain) in enumerate(all_strains)
                strain in ctx.undone_strains && continue
                i > length(record.sample_data) && continue
                fmt = parse_format_field(record.format_keys, record.sample_data[i])
                dp_str = get(fmt, "DP", "0")
                dp = isempty(dp_str) || dp_str == "." ? 0 : parse(Int, dp_str)
                if dp > 0
                    prev_coverage_span[strain] = (record.chrom, record.end_pos, dp)
                end
            end

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
                (_, _, ref, alt, cann_str) = parsed
                cache_entries[(ref, alt)] = CacheEntry(cann_str)
            end
            advance!(cache_pf)
        end

        if handle_variant_record!(record, cache_entries, ctx, writers, transcript_cache, all_strains, prev_coverage_span)
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
