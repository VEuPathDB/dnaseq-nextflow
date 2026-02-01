#!/usr/bin/env julia

# processSequenceVariations.jl
# Replaces bin/processSequenceVariationsNew.pl
#
# Reads a sorted SNP file and a sorted cache file, merges them in a single
# pass, annotates coding variants with codon/product information via SQLite
# sequence and indel databases, and writes four output files:
#   cache (variationFeature.dat after rename), snpFeature.dat, allele.dat, product.dat

using SQLite
using Printf

# ---------------------------------------------------------------------------
# CLI argument parsing (hand-rolled, no package dependency)
# ---------------------------------------------------------------------------

function parse_args(args::Vector{String})
    flags = Dict{String,String}()
    i = 1
    while i <= length(args)
        if startswith(args[i], "--")
            key = args[i][3:end]  # strip leading --
            if i + 1 <= length(args) && !startswith(args[i+1], "--")
                flags[key] = args[i+1]
                i += 2
            else
                flags[key] = ""   # boolean flag
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
# IUPAC ambiguity expansion
# ---------------------------------------------------------------------------

const IUPAC_EXPAND = Dict{Char,Vector{Char}}(
    'A'=>['A'], 'G'=>['G'], 'C'=>['C'], 'T'=>['T'],
    'R'=>['A','G'], 'Y'=>['C','T'], 'K'=>['G','T'], 'M'=>['A','C'],
    'S'=>['G','C'], 'W'=>['A','T'],
    'B'=>['G','T','C'], 'D'=>['G','A','T'], 'H'=>['A','C','T'], 'V'=>['G','C','A'],
    'N'=>['A','G','C','T']
)

"""
    expand_codon(codon::String) -> Vector{String}

Cartesian product expansion of IUPAC ambiguity codes in a 3-character codon.
If the codon is not exactly 3 characters, returns ["NNN"] expanded (i.e. all 64 codons).
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

"""
    translate_codon(codon::String) -> String
"""
function translate_codon(codon::String)
    get(CODON_TABLE, uppercase(codon), "X")
end

# ---------------------------------------------------------------------------
# GTF parsing → CDS interval structures
# ---------------------------------------------------------------------------

"""
    struct CDSInterval
        seq_id::String          # chromosome / sequence source id
        start_pos::Int          # 1-based inclusive start (genomic)
        end_pos::Int            # 1-based inclusive end (genomic)
        strand::Int             # +1 or -1
        transcript_id::String
        cds_number::Int         # exon order within this transcript (1-based)
    end
"""
struct CDSInterval
    seq_id::String
    start_pos::Int
    end_pos::Int
    strand::Int
    transcript_id::String
    cds_number::Int
end

"""
    TranscriptInfo holds the ordered CDS exons for one transcript.
    exons are stored in genomic order (sorted by start_pos).
"""
struct TranscriptInfo
    seq_id::String
    strand::Int
    exons::Vector{CDSInterval}   # sorted by start_pos
end

"""
    parse_gtf(gtf_file::String) -> (Vector{CDSInterval}, Dict{String, TranscriptInfo})

Parse a GTF file, extracting CDS features. Returns:
  - sorted array of all CDS intervals (for binary search)
  - per-transcript info dict
"""
function parse_gtf(gtf_file::String)
    # First pass: collect CDS entries grouped by transcript
    transcript_cds = Dict{String, Vector{Tuple{String,Int,Int,Int}}}()  # tid -> [(seq_id, start, end, strand)]
    transcript_order = String[]  # preserve first-seen order

    open(gtf_file, "r") do fh
        for line in eachline(fh)
            startswith(line, '#') && continue
            fields = split(line, '\t')
            length(fields) < 9 && continue
            feature = fields[3]
            feature == "CDS" || continue

            seq_id = fields[1]
            start_pos = parse(Int, fields[4])
            end_pos = parse(Int, fields[5])
            strand_str = fields[7]
            strand = (strand_str == "-") ? -1 : 1
            attributes = fields[9]

            # Extract transcript_id from attributes
            tid = ""
            for attr in split(attributes, ';')
                attr = strip(attr)
                if startswith(attr, "transcript_id")
                    # transcript_id "VALUE" or transcript_id VALUE
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

    # Build structures
    all_intervals = CDSInterval[]
    transcript_info = Dict{String, TranscriptInfo}()

    for tid in transcript_order
        entries = transcript_cds[tid]
        # Sort exons by start position (genomic order)
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

    # Sort all_intervals by (seq_id, start_pos) for binary search
    sort!(all_intervals; by = iv -> (iv.seq_id, iv.start_pos))

    (all_intervals, transcript_info)
end

# ---------------------------------------------------------------------------
# Binary search: find CDS interval containing a genomic location
# Returns nothing if not in any CDS
# ---------------------------------------------------------------------------

function find_cds(intervals::Vector{CDSInterval}, seq_id::String, location::Int)
    # Binary search for the rightmost interval with start_pos <= location on matching seq_id
    lo, hi = 1, length(intervals)
    candidate = 0   # index of rightmost interval with (seq_id, start_pos) <= (seq_id, location)

    while lo <= hi
        mid = (lo + hi) ÷ 2
        iv = intervals[mid]
        # Compare (seq_id, start_pos) vs (seq_id, location)
        if iv.seq_id < seq_id || (iv.seq_id == seq_id && iv.start_pos <= location)
            candidate = mid
            lo = mid + 1
        else
            hi = mid - 1
        end
    end

    if candidate == 0
        return nothing
    end

    iv = intervals[candidate]
    if iv.seq_id == seq_id && location >= iv.start_pos && location <= iv.end_pos
        return iv
    end

    nothing
end

# ---------------------------------------------------------------------------
# Compute position_in_cds from transcript info + genomic location
# position_in_cds is 1-based: sum of lengths of all prior exons + offset in current exon
# ---------------------------------------------------------------------------

function compute_position_in_cds(tinfo::TranscriptInfo, location::Int, cds_number::Int)
    prior_len = 0
    for (idx, exon) in enumerate(tinfo.exons)
        if idx == cds_number
            # offset within this exon (1-based)
            return prior_len + (location - exon.start_pos) + 1
        end
        # exon length = end - start + 1, but Perl uses end - start (no +1 for prior lengths)
        # Matching Perl: prior_ref_cds_len += cds_end - cds_start
        prior_len += exon.end_pos - exon.start_pos
    end
    # Should not reach here if cds_number is valid
    error("cds_number $cds_number out of range for transcript")
end

# ---------------------------------------------------------------------------
# position_in_codon (called position_in_protein in product file)
# Matches Perl: (pos % 3 == 0) ? 3 : (pos % 3)
# ---------------------------------------------------------------------------

function position_in_codon(pos_in_cds::Int)
    r = pos_in_cds % 3
    r == 0 ? 3 : r
end

# ---------------------------------------------------------------------------
# Extract codon from a sequence given a 1-based position_in_cds
# ---------------------------------------------------------------------------

function extract_codon(sequence::String, pos_in_cds::Int)
    pic = position_in_codon(pos_in_cds)
    # codon_start is 0-based index of first base of the codon
    codon_start = pos_in_cds - pic  # 0-based
    if codon_start + 3 > length(sequence)
        return "NNN"   # incomplete codon at end of CDS
    end
    sequence[codon_start+1 : codon_start+3]
end

# ---------------------------------------------------------------------------
# Frameshift precomputation
# ---------------------------------------------------------------------------

"""
    frameshift_info::Dict{String, Dict{String, Tuple{Bool, Int}}}
    strain -> transcript -> (has_frameshift, frameshift_position)

    frameshift_position is the genomic position of the first indel in that
    transcript where the cumulative shift % 3 != 0.
"""

function precompute_frameshifts(indel_db::SQLite.DB, transcript_info::Dict{String,TranscriptInfo})
    # For each (strain, transcript) that has indels, compute running cumulative sum
    # and record first position where cumsum % 3 != 0
    fs_info = Dict{String, Dict{String, Tuple{Bool, Int}}}()

    # Get all distinct (strain, transcript_id) pairs that have indels
    # For each, get indels sorted by position
    stmt = SQLite.Stmt(indel_db,
        "SELECT strain, transcript_id, position, shift_amount FROM indels ORDER BY strain, transcript_id, position")

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

    for row in SQLite.rows(stmt)
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

    fs_info
end

# ---------------------------------------------------------------------------
# Check downstream-of-frameshift for a variant
# ---------------------------------------------------------------------------

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

    # Per plan: transcript coordinates are always 5'→3', so comparison is always >
    location > fs_pos ? 1 : 0
end

# ---------------------------------------------------------------------------
# SQLite sequence loading
# ---------------------------------------------------------------------------

"""
    load_transcript_sequences(db, transcript_id) -> Dict{String,String}
    Fetches all strain sequences for a transcript in one query.
"""
function load_transcript_sequences(db::SQLite.DB, transcript_id::String)
    stmt = SQLite.Stmt(db, "SELECT strain, sequence FROM sequences WHERE transcript_id = ?")
    result = Dict{String,String}()
    for row in SQLite.rows(stmt, transcript_id)
        result[row[1]::String] = row[2]::String
    end
    result
end

"""
    get_indel_shift(db, transcript_id, strain, position) -> Int
    Returns sum of shift_amount for indels before the given position.
"""
function get_indel_shift(db::SQLite.DB, transcript_id::String, strain::String, position::Int)
    stmt = SQLite.Stmt(db,
        "SELECT COALESCE(SUM(shift_amount), 0) FROM indels WHERE transcript_id = ? AND strain = ? AND position < ?")
    row = first(SQLite.rows(stmt, transcript_id, strain, position))
    row[1]::Int
end

# ---------------------------------------------------------------------------
# Coverage file reader
# ---------------------------------------------------------------------------

"""
    CoverageReader wraps an open coverage file with a peek buffer.
    Coverage files are tab-separated: seq_id, start, end, coverage_array, percent_array
"""
struct CoverageReader
    fh::IOStream
    # Mutable state via a 1-element vector (workaround for immutable struct)
    state::Vector{Any}  # [peek_fields, coverage_array, percents_array, exhausted]
end

function open_coverage_reader(path::String)
    fh = open(path, "r")
    cr = CoverageReader(fh, Any[nothing, nothing, nothing, false])
    # Read first line into peek
    advance_coverage!(cr)
    cr
end

function advance_coverage!(cr::CoverageReader)
    line = readline(cr.fh)
    if isempty(line) && eof(cr.fh)
        cr.state[1] = nothing
        cr.state[4] = true   # exhausted
    else
        cr.state[1] = split(line, '\t')
        cr.state[2] = nothing   # clear cached arrays
        cr.state[3] = nothing
    end
end

function coverage_exhausted(cr::CoverageReader)
    cr.state[4] == true
end

function coverage_peek(cr::CoverageReader)
    cr.state[1]  # returns the split fields or nothing
end

"""
    get_coverage(cr, seq_id, location) -> (coverage, percent) or nothing
    Advances the reader as needed. Returns nothing if position not covered.
"""
function get_coverage(cr::CoverageReader, seq_id::String, location::Int)
    while !coverage_exhausted(cr)
        peek = coverage_peek(cr)
        isnothing(peek) && break

        p_seq_id = peek[1]
        p_start = parse(Int, peek[2])
        p_end = parse(Int, peek[3])

        # If current span covers the location, extract values
        if p_seq_id == seq_id && location >= p_start && location <= p_end
            # Lazily parse coverage/percent arrays
            if isnothing(cr.state[2])
                cr.state[2] = split(peek[4], ',')   # coverage array (strings)
                cr.state[3] = split(peek[5], ',')   # percent array (strings)
            end
            idx = location - p_start + 1  # 1-based index
            cov = cr.state[2][idx]
            pct = cr.state[3][idx]
            return (cov, pct)
        end

        # If we've passed the location, stop
        if p_seq_id > seq_id || (p_seq_id == seq_id && p_start > location)
            break
        end

        # Advance to next span
        advance_coverage!(cr)
    end
    nothing
end

function close_coverage_reader(cr::CoverageReader)
    close(cr.fh)
end

# ---------------------------------------------------------------------------
# Open all coverage files in a directory
# ---------------------------------------------------------------------------

function open_coverage_files(coverage_dir::String)
    readers = Dict{String, CoverageReader}()
    for fname in readdir(coverage_dir)
        if endswith(fname, ".coverage.txt")
            strain = replace(fname, r"\.coverage\.txt$" => "")
            path = joinpath(coverage_dir, fname)
            readers[strain] = open_coverage_reader(path)
        end
    end
    readers
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
    product::Vector{String}     # list of translated products
    reference_codon::String
    codon::String
    has_nonsynonomous::Int      # typo preserved per spec
    cds_number::Int
    matches_reference::Int
end

function Variation()
    Variation("", 0, "", "", "", "", "", "", "", "",
              0, 0, 0, 0, "", String[], "", "", 0, 0, 0)
end

# ---------------------------------------------------------------------------
# Parse a cache line (20 columns in new format) into a Variation
# Cache columns: sequence_source_id, location, strain, reference, base,
#   coverage, percent, quality, pvalue, snp_source_id, is_coding,
#   position_in_cds, position_in_codon, downstream_of_frameshift,
#   transcript, product, reference_codon, codon, has_nonsynonomous, cds_number
# ---------------------------------------------------------------------------

function parse_cache_line(line::String)
    fields = split(line, '\t')
    length(fields) < 20 && return nothing

    v = Variation()
    v.sequence_source_id = fields[1]
    v.location = parse(Int, fields[2])
    v.strain = fields[3]
    v.reference = fields[4]
    v.base = fields[5]
    v.coverage = fields[6]
    v.percent = fields[7]
    v.quality = fields[8]
    v.pvalue = fields[9]
    v.snp_source_id = fields[10]
    v.is_coding = parse(Int, fields[11])
    v.position_in_cds = parse(Int, fields[12])
    v.position_in_codon = parse(Int, fields[13])
    v.downstream_of_frameshift = parse(Int, fields[14])
    v.transcript = fields[15]
    v.product = split(fields[16], ':')
    v.reference_codon = fields[17]
    v.codon = fields[18]
    v.has_nonsynonomous = parse(Int, fields[19])
    v.cds_number = parse(Int, fields[20])
    v
end

# ---------------------------------------------------------------------------
# Parse a SNP file line (10 columns) into a Variation
# SNP columns: seqId, location, strain, ref, base, coverage, percent, quality, pvalue, snp_source_id
# ---------------------------------------------------------------------------

function parse_snp_line(line::String)
    fields = split(line, '\t')
    length(fields) < 10 && return nothing

    v = Variation()
    v.sequence_source_id = fields[1]
    v.location = parse(Int, fields[2])
    v.strain = fields[3]
    v.reference = fields[4]
    v.base = fields[5]
    v.coverage = fields[6]
    v.percent = fields[7]
    v.quality = fields[8]
    v.pvalue = fields[9]
    v.snp_source_id = fields[10]
    v
end

# ---------------------------------------------------------------------------
# Sorted-merge iteration helpers
# ---------------------------------------------------------------------------

"""
    PeekedFile wraps a file handle with a one-line lookahead.
    line: current peeked line (empty string if exhausted)
    exhausted: true when EOF reached
"""
mutable struct PeekedFile
    fh::IOStream
    line::String
    exhausted::Bool
end

function open_peeked(path::String)
    fh = open(path, "r")
    pf = PeekedFile(fh, "", false)
    advance!(pf)
    pf
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

# Extract (seq_id, location) from a peeked line for sorting
# Cache line: col 0 = seq_id, col 1 = location
# SNP line: col 0 = seq_id, col 1 = location
function peek_sort_key(line::String)
    # Find first two tab-separated fields
    tab1 = findfirst('\t', line)
    isnothing(tab1) && return ("", 0)
    seq_id = line[1:tab1-1]
    rest = line[tab1+1:end]
    tab2 = findfirst('\t', rest)
    loc_str = isnothing(tab2) ? rest : rest[1:tab2-1]
    (seq_id, parse(Int, loc_str))
end

# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------

function main()
    args = parse_args(ARGS)

    snp_file = args["snp_file"]
    cache_file = args["cache_file"]
    undone_strains_file = args["undone_strains_file"]
    coverage_dir = args["coverage_directory"]
    transcript_db_path = args["transcript_db"]
    indel_db_path = args["indel_db"]
    gtf_file = args["gtf_file"]
    reference_strain = args["reference_strain"]

    # Read undone strains
    undone_strains = Set{String}()
    if isfile(undone_strains_file)
        open(undone_strains_file, "r") do fh
            for line in eachline(fh)
                s = strip(line)
                !isempty(s) && push!(undone_strains, s)
            end
        end
    end

    # Parse GTF
    (cds_intervals, transcript_info) = parse_gtf(gtf_file)

    # Open SQLite databases (read-only)
    transcript_db = SQLite.DB(transcript_db_path; readonly=true)
    indel_db = SQLite.DB(indel_db_path; readonly=true)

    # Precompute frameshift info
    fs_info = precompute_frameshifts(indel_db, transcript_info)

    # Open coverage files
    coverage_readers = open_coverage_files(coverage_dir)
    all_strains = collect(keys(coverage_readers))

    # Determine output directory (same as cache file's directory)
    cache_dir = dirname(abspath(cache_file))
    temp_cache_file = joinpath(cache_dir, "cache.tmp")
    snp_output_file = joinpath(cache_dir, "snpFeature.dat")
    allele_output_file = joinpath(cache_dir, "allele.dat")
    product_output_file = joinpath(cache_dir, "product.dat")

    # Open output files
    cache_fh = open(temp_cache_file, "w")
    snp_fh = open(snp_output_file, "w")
    allele_fh = open(allele_output_file, "w")
    product_fh = open(product_output_file, "w")

    # Open input files with peek
    cache_pf = open_peeked(cache_file)
    snp_pf = open_peeked(snp_file)

    # Sequence cache: cleared when transcript changes
    current_transcript_id = ""
    current_transcript_seqs = Dict{String,String}()

    count = 0

    # ---------------------------------------------------------------------------
    # Main loop: sorted merge
    # ---------------------------------------------------------------------------
    while !cache_pf.exhausted || !snp_pf.exhausted
        # Determine the next (seq_id, location) to process
        cache_key = cache_pf.exhausted ? ("~", typemax(Int)) : peek_sort_key(cache_pf.line)
        snp_key = snp_pf.exhausted ? ("~", typemax(Int)) : peek_sort_key(snp_pf.line)

        # Pick the smaller key (seq_id first, then location)
        if cache_key[1] < snp_key[1] || (cache_key[1] == snp_key[1] && cache_key[2] <= snp_key[2])
            current_key = cache_key
        else
            current_key = snp_key
        end

        current_seq_id, current_location = current_key

        # Collect all variations at this position
        variations = Variation[]
        cache_strains = Set{String}()

        # Drain cache lines at this position
        while !cache_pf.exhausted
            k = peek_sort_key(cache_pf.line)
            k[1] != current_seq_id && break
            k[2] != current_location && break

            v = parse_cache_line(cache_pf.line)
            advance!(cache_pf)
            isnothing(v) && continue

            # Skip undone strains
            if v.strain in undone_strains
                continue
            end

            push!(variations, v)
            push!(cache_strains, v.strain)
        end

        # Drain SNP lines at this position (skip strains already from cache)
        while !snp_pf.exhausted
            k = peek_sort_key(snp_pf.line)
            k[1] != current_seq_id && break
            k[2] != current_location && break

            v = parse_snp_line(snp_pf.line)
            advance!(snp_pf)
            isnothing(v) && continue

            # Skip if strain already present from cache
            if v.strain in cache_strains
                continue
            end

            push!(variations, v)
        end

        # Skip if empty batch
        isempty(variations) && continue

        # ---------------------------------------------------------------------------
        # Annotate: CDS lookup
        # ---------------------------------------------------------------------------
        cds_hit = find_cds(cds_intervals, current_seq_id, current_location)

        is_coding = 0
        transcript_id = ""
        cds_number = 0
        pos_in_cds = 0
        pos_in_codon_val = 0

        if !isnothing(cds_hit)
            is_coding = 1
            transcript_id = cds_hit.transcript_id
            cds_number = cds_hit.cds_number

            tinfo = transcript_info[transcript_id]
            pos_in_cds = compute_position_in_cds(tinfo, current_location, cds_number)
            pos_in_codon_val = position_in_codon(pos_in_cds)
        end

        # Load transcript sequences if transcript changed
        if is_coding == 1 && transcript_id != current_transcript_id
            current_transcript_id = transcript_id
            current_transcript_seqs = load_transcript_sequences(transcript_db, transcript_id)
        end

        # ---------------------------------------------------------------------------
        # Compute reference codon/product
        # ---------------------------------------------------------------------------
        ref_codon = ""
        ref_product = ""

        if is_coding == 1
            ref_seq = get(current_transcript_seqs, reference_strain, "")
            if !isempty(ref_seq)
                ref_codon = extract_codon(ref_seq, pos_in_cds)
                ref_product = translate_codon(ref_codon)
            end
        end

        # ---------------------------------------------------------------------------
        # Annotate each variation with coding info
        # ---------------------------------------------------------------------------
        for v in variations
            # Skip variations that already have position_in_codon set (from cache)
            if v.position_in_codon != 0
                # Already annotated (came from cache) — keep existing values
                # but update is_coding/transcript if needed
                continue
            end

            v.is_coding = is_coding
            v.transcript = transcript_id
            v.cds_number = cds_number

            if is_coding == 1
                v.position_in_cds = pos_in_cds
                v.position_in_codon = pos_in_codon_val
                v.reference_codon = ref_codon

                # Get strain's sequence and compute adjusted position
                strain_seq = get(current_transcript_seqs, v.strain, "")
                if !isempty(strain_seq)
                    # Compute position adjustment via indel DB
                    shift = get_indel_shift(indel_db, transcript_id, v.strain, pos_in_cds)
                    adjusted_pos = pos_in_cds + shift

                    strain_codon = extract_codon(strain_seq, adjusted_pos)
                    v.codon = strain_codon

                    # Expand IUPAC and translate
                    expanded = expand_codon(strain_codon)
                    v.product = [translate_codon(c) for c in expanded]
                else
                    v.codon = "NNN"
                    v.product = ["X"]
                end
            end
        end

        # ---------------------------------------------------------------------------
        # Build reference variation record
        # ---------------------------------------------------------------------------
        ref_allele = variations[1].reference
        # If reference field is empty (SNP-only lines), use the base from the first entry
        if isempty(ref_allele)
            ref_allele = variations[1].base
        end

        adjacent_snp_causes_product_difference = 0
        if is_coding == 1
            for v in variations
                product_str = join(v.product, ":")
                if product_str != ref_product
                    adjacent_snp_causes_product_difference = 1
                    break
                end
            end
        end

        # Mark has_nonsynonomous on each variation
        for v in variations
            if is_coding == 1
                product_str = join(v.product, ":")
                v.has_nonsynonomous = (product_str != ref_product) ? 1 : 0
            end
        end

        ref_variation = Variation(
            current_seq_id,
            current_location,
            reference_strain,
            ref_allele,
            ref_allele,       # base = reference allele
            "", "",           # coverage, percent (empty for reference)
            "", "",           # quality, pvalue
            "",               # snp_source_id
            is_coding,
            pos_in_cds,
            pos_in_codon_val,
            0,                # downstream_of_frameshift (computed later)
            transcript_id,
            is_coding == 1 ? [ref_product] : String[],
            ref_codon,
            ref_codon,        # codon = ref_codon for reference
            adjacent_snp_causes_product_difference,
            cds_number,
            1                 # matches_reference
        )

        # Add reference variation to the batch
        push!(variations, ref_variation)

        # ---------------------------------------------------------------------------
        # Check if there is actual variation (more than one distinct allele)
        # ---------------------------------------------------------------------------
        alleles_set = Set(v.base for v in variations)
        if length(alleles_set) <= 1
            continue  # no variation, skip
        end

        # ---------------------------------------------------------------------------
        # Fill coverage for strains not in the variation batch
        # ---------------------------------------------------------------------------
        variation_strains = Set(v.strain for v in variations)

        for strain in all_strains
            strain in variation_strains && continue

            cr = coverage_readers[strain]
            cov_result = get_coverage(cr, current_seq_id, current_location)
            if !isnothing(cov_result)
                cv = Variation()
                cv.sequence_source_id = current_seq_id
                cv.location = current_location
                cv.strain = strain
                cv.reference = ref_allele
                cv.base = ref_allele
                cv.coverage = cov_result[1]
                cv.percent = cov_result[2]
                cv.matches_reference = 1
                cv.is_coding = is_coding
                cv.transcript = transcript_id
                cv.position_in_cds = pos_in_cds
                cv.position_in_codon = pos_in_codon_val
                cv.reference_codon = ref_codon
                cv.codon = ref_codon
                cv.product = is_coding == 1 ? [ref_product] : String[]
                cv.has_nonsynonomous = adjacent_snp_causes_product_difference
                cv.cds_number = cds_number
                push!(variations, cv)
            end
        end

        # ---------------------------------------------------------------------------
        # Compute downstream_of_frameshift and write cache
        # ---------------------------------------------------------------------------
        for v in variations
            if !isempty(v.transcript)
                v.downstream_of_frameshift = check_downstream_of_frameshift(
                    fs_info, v.strain, v.transcript, current_location)
            else
                v.downstream_of_frameshift = 0
            end

            # Skip reference strain in cache output
            v.strain == reference_strain && continue

            # matches_reference
            v.matches_reference = (v.base == ref_allele) ? 1 : 0

            # Write cache line (20 columns)
            # sequence_source_id, location, strain, reference, base, coverage, percent,
            # quality, pvalue, snp_source_id, is_coding, position_in_cds, position_in_codon,
            # downstream_of_frameshift, transcript, product, reference_codon, codon,
            # has_nonsynonomous, cds_number
            product_str = join(v.product, ":")
            write(cache_fh, join([
                v.sequence_source_id,
                string(v.location),
                v.strain,
                v.reference,
                v.base,
                v.coverage,
                v.percent,
                v.quality,
                v.pvalue,
                v.snp_source_id,
                string(v.is_coding),
                string(v.position_in_cds),
                string(v.position_in_codon),
                string(v.downstream_of_frameshift),
                v.transcript,
                product_str,
                v.reference_codon,
                v.codon,
                string(v.has_nonsynonomous),
                string(v.cds_number)
            ], "\t"), "\n")
        end

        # ---------------------------------------------------------------------------
        # Write SNP feature summary
        # ---------------------------------------------------------------------------
        # Columns: location, transcript_id, source_id, reference_strain, reference_na,
        #   has_nonsynonymous_allele, major_allele, minor_allele, major_allele_count,
        #   minor_allele_count, major_product, minor_product, distinct_strain_count,
        #   distinct_allele_count, has_coding_mutation, total_allele_count,
        #   has_stop_codon, ref_codon

        allele_counts = Dict{String,Int}()
        product_counts = Dict{String,Int}()
        strain_set = Set{String}()
        has_stop_codon = 0
        total_allele_count = length(variations)

        for v in variations
            allele_counts[v.base] = get(allele_counts, v.base, 0) + 1
            strain_set |= Set([v.strain])

            # Products are stored as colon-joined string in cache output;
            # for SNP feature we split and count individual products
            for p in v.product
                product_counts[p] = get(product_counts, p, 0) + 1
                p == "*" && (has_stop_codon = 1)
            end
        end

        distinct_strain_count = length(strain_set)
        distinct_allele_count = length(allele_counts)
        has_nonsynonymous_allele = length(product_counts) > 1 ? 1 : 0

        # Sort alleles: descending by count, then ascending alphabetically
        sorted_alleles = sort(collect(keys(allele_counts));
            by = a -> (-allele_counts[a], a))
        sorted_products = sort(collect(keys(product_counts));
            by = p -> (-product_counts[p], p))

        major_allele = sorted_alleles[1]
        minor_allele = length(sorted_alleles) > 1 ? sorted_alleles[2] : ""
        major_allele_count = allele_counts[major_allele]
        minor_allele_count = length(sorted_alleles) > 1 ? allele_counts[minor_allele] : ""

        major_product = length(sorted_products) > 0 ? sorted_products[1] : ""
        minor_product = length(sorted_products) > 1 ? sorted_products[2] : ""

        write(snp_fh, join([
            string(current_location),
            transcript_id,
            current_seq_id,
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
            string(is_coding),
            string(total_allele_count),
            string(has_stop_codon),
            ref_codon
        ], "\t"), "\n")

        # ---------------------------------------------------------------------------
        # Write allele and product files (only for coding variants)
        # ---------------------------------------------------------------------------
        if is_coding == 1
            # Allele file: allele, distinct_strain_count, allele_count, average_coverage, average_read_percent
            for allele in keys(allele_counts)
                distinct_strains_for_allele = Set{String}()
                allele_count = 0
                sum_coverage = 0.0
                sum_percent = 0.0

                for v in variations
                    v.base != allele && continue
                    allele_count += 1
                    push!(distinct_strains_for_allele, v.strain)
                    # coverage and percent may be empty strings
                    cov = isempty(v.coverage) ? 0.0 : parse(Float64, v.coverage)
                    pct = isempty(v.percent) ? 0.0 : parse(Float64, v.percent)
                    sum_coverage += cov
                    sum_percent += pct
                end

                avg_cov = allele_count > 0 ? sum_coverage / allele_count : 0.0
                avg_pct = allele_count > 0 ? sum_percent / allele_count : 0.0

                write(allele_fh, join([
                    allele,
                    string(length(distinct_strains_for_allele)),
                    string(allele_count),
                    @sprintf("%.2f", avg_cov),
                    @sprintf("%.2f", avg_pct)
                ], "\t"), "\n")
            end

            # Product file: codon, position_in_protein, transcript, count, product,
            #   ref_location_cds, ref_location_protein, downstream_of_frameshift
            #
            # One row per expanded (non-ambiguous) codon per variant.
            # Recount products across all variations for the count column.
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
                    count = get(all_product_counts, product, 0)

                    write(product_fh, join([
                        ec,
                        string(pos_in_codon_val),
                        transcript_id,
                        string(count),
                        product,
                        string(pos_in_cds),
                        string(pos_in_codon_val),
                        string(v.downstream_of_frameshift)
                    ], "\t"), "\n")
                end
            end
        end

        count += 1
        if count % 1000 == 0
            @info "Processed $count SNPs"
        end
    end

    # ---------------------------------------------------------------------------
    # Close files and finalize
    # ---------------------------------------------------------------------------
    close(cache_fh)
    close(snp_fh)
    close(allele_fh)
    close(product_fh)
    close_peeked(cache_pf)
    close_peeked(snp_pf)

    for cr in values(coverage_readers)
        close_coverage_reader(cr)
    end

    close(transcript_db)
    close(indel_db)

    # Rename temp cache to final cache
    if isfile(cache_file)
        rm(cache_file)
    end
    mv(temp_cache_file, cache_file)

    # Truncate the SNP input file (consumed)
    open(snp_file, "w") do fh end

    # Truncate the undone strains file (consumed)
    open(undone_strains_file, "w") do fh end
end

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
main()
