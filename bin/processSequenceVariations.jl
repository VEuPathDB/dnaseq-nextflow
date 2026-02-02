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
# Processing data structures for refactored main()
# ---------------------------------------------------------------------------

"""
    ProcessingContext holds all read-only reference data and database connections.
"""
struct ProcessingContext
    reference_strain::String
    undone_strains::Set{String}
    cds_intervals::Vector{CDSInterval}
    transcript_info::Dict{String,TranscriptInfo}
    transcript_db::SQLite.DB
    indel_db::SQLite.DB
    fs_info::Dict{String, Dict{String, Tuple{Bool, Int}}}  # frameshift info
    coverage_readers::Dict{String, CoverageReader}
    all_strains::Vector{String}
end

"""
    OutputWriters encapsulates the four output file handles.
"""
struct OutputWriters
    cache_fh::IOStream
    snp_fh::IOStream
    allele_fh::IOStream
    product_fh::IOStream
end

"""
    PositionAnnotation bundles all annotation data computed for a genomic position.
"""
struct PositionAnnotation
    is_coding::Int
    transcript_id::String
    cds_number::Int
    pos_in_cds::Int
    pos_in_codon_val::Int
    ref_codon::String
    ref_product::String
end

"""
    TranscriptSequenceCache explicitly manages mutable state for transcript sequence caching.
"""
mutable struct TranscriptSequenceCache
    current_transcript_id::String
    current_transcript_seqs::Dict{String,String}
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
# Resource management functions
# ---------------------------------------------------------------------------

"""
    initialize_processing_context(args) -> ProcessingContext

Initialize all read-only reference data and database connections.
Loads undone strains, parses GTF, opens databases, precomputes frameshifts, opens coverage files.
"""
function initialize_processing_context(args)
    # Read undone strains
    undone_strains = Set{String}()
    undone_strains_file = args["undone_strains_file"]
    if isfile(undone_strains_file)
        open(undone_strains_file, "r") do fh
            for line in eachline(fh)
                s = strip(line)
                !isempty(s) && push!(undone_strains, s)
            end
        end
    end

    # Parse GTF
    (cds_intervals, transcript_info) = parse_gtf(args["gtf_file"])

    # Open SQLite databases (read-only)
    transcript_db = SQLite.DB(args["transcript_db"]; readonly=true)
    indel_db = SQLite.DB(args["indel_db"]; readonly=true)

    # Precompute frameshift info
    fs_info = precompute_frameshifts(indel_db, transcript_info)

    # Open coverage files
    coverage_readers = open_coverage_files(args["coverage_directory"])
    all_strains = collect(keys(coverage_readers))

    ProcessingContext(
        args["reference_strain"],
        undone_strains,
        cds_intervals,
        transcript_info,
        transcript_db,
        indel_db,
        fs_info,
        coverage_readers,
        all_strains
    )
end

"""
    open_output_writers(cache_file) -> (OutputWriters, temp_cache_file)

Creates output file handles and returns temp cache path.
"""
function open_output_writers(cache_file)
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

    (OutputWriters(cache_fh, snp_fh, allele_fh, product_fh), temp_cache_file)
end

"""
    close_processing_context(ctx::ProcessingContext)

Closes coverage readers and databases.
"""
function close_processing_context(ctx::ProcessingContext)
    for cr in values(ctx.coverage_readers)
        close_coverage_reader(cr)
    end

    close(ctx.transcript_db)
    close(ctx.indel_db)
end

"""
    close_output_writers(writers::OutputWriters)

Closes all output file handles.
"""
function close_output_writers(writers::OutputWriters)
    close(writers.cache_fh)
    close(writers.snp_fh)
    close(writers.allele_fh)
    close(writers.product_fh)
end

"""
    finalize_output_files(cache_file, temp_cache_file, snp_file, undone_strains_file)

Renames temp cache, truncates consumed input files.
"""
function finalize_output_files(cache_file, temp_cache_file, snp_file, undone_strains_file)
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
# Position processing functions
# ---------------------------------------------------------------------------

"""
    determine_next_position(cache_pf, snp_pf) -> (seq_id, location) or nothing

Determines which position to process next from peeked files.
Returns nothing if both files are exhausted.
"""
function determine_next_position(cache_pf::PeekedFile, snp_pf::PeekedFile)
    if cache_pf.exhausted && snp_pf.exhausted
        return nothing
    end

    # Determine the next (seq_id, location) to process
    cache_key = cache_pf.exhausted ? ("~", typemax(Int)) : peek_sort_key(cache_pf.line)
    snp_key = snp_pf.exhausted ? ("~", typemax(Int)) : peek_sort_key(snp_pf.line)

    # Pick the smaller key (seq_id first, then location)
    if cache_key[1] < snp_key[1] || (cache_key[1] == snp_key[1] && cache_key[2] <= snp_key[2])
        cache_key
    else
        snp_key
    end
end

"""
    collect_variations_at_position(cache_pf, snp_pf, seq_id, location, undone_strains)
        -> (variations, cache_strains)

Drains cache and SNP files for current position.
Filters undone strains.
"""
function collect_variations_at_position(
    cache_pf::PeekedFile,
    snp_pf::PeekedFile,
    current_seq_id::String,
    current_location::Int,
    undone_strains::Set{String}
)
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

    (variations, cache_strains)
end

"""
    has_variation(variations) -> bool

Checks if there's actual variation (>1 distinct allele).
"""
function has_variation(variations::Vector{Variation})
    alleles_set = Set(v.base for v in variations)
    length(alleles_set) > 1
end

# ---------------------------------------------------------------------------
# Annotation logic functions
# ---------------------------------------------------------------------------

"""
    annotate_position(seq_id, location, ctx, transcript_cache) -> PositionAnnotation

CDS lookup, position in CDS/codon computation.
Loads transcript sequences, computes reference codon/product.
Updates transcript_cache if transcript changes.
"""
function annotate_position(
    seq_id::String,
    location::Int,
    ctx::ProcessingContext,
    transcript_cache::TranscriptSequenceCache
)
    # CDS lookup
    cds_hit = find_cds(ctx.cds_intervals, seq_id, location)

    is_coding = 0
    transcript_id = ""
    cds_number = 0
    pos_in_cds = 0
    pos_in_codon_val = 0
    ref_codon = ""
    ref_product = ""

    if !isnothing(cds_hit)
        is_coding = 1
        transcript_id = cds_hit.transcript_id
        cds_number = cds_hit.cds_number

        tinfo = ctx.transcript_info[transcript_id]
        pos_in_cds = compute_position_in_cds(tinfo, location, cds_number)
        pos_in_codon_val = position_in_codon(pos_in_cds)

        # Load transcript sequences if transcript changed
        if transcript_id != transcript_cache.current_transcript_id
            transcript_cache.current_transcript_id = transcript_id
            transcript_cache.current_transcript_seqs = load_transcript_sequences(ctx.transcript_db, transcript_id)
        end

        # Compute reference codon/product
        ref_seq = get(transcript_cache.current_transcript_seqs, ctx.reference_strain, "")
        if !isempty(ref_seq)
            ref_codon = extract_codon(ref_seq, pos_in_cds)
            ref_product = translate_codon(ref_codon)
        end
    end

    PositionAnnotation(
        is_coding,
        transcript_id,
        cds_number,
        pos_in_cds,
        pos_in_codon_val,
        ref_codon,
        ref_product
    )
end

"""
    annotate_variations!(variations, annotation, ctx, transcript_cache)

Annotates each variation with coding info.
Mutates variation records in place.
"""
function annotate_variations!(
    variations::Vector{Variation},
    annotation::PositionAnnotation,
    ctx::ProcessingContext,
    transcript_cache::TranscriptSequenceCache
)
    for v in variations
        # Skip variations that already have position_in_codon set (from cache)
        if v.position_in_codon != 0
            # Already annotated (came from cache) — keep existing values
            continue
        end

        v.is_coding = annotation.is_coding
        v.transcript = annotation.transcript_id
        v.cds_number = annotation.cds_number

        if annotation.is_coding == 1
            v.position_in_cds = annotation.pos_in_cds
            v.position_in_codon = annotation.pos_in_codon_val
            v.reference_codon = annotation.ref_codon

            # Get strain's sequence and compute adjusted position
            strain_seq = get(transcript_cache.current_transcript_seqs, v.strain, "")
            if !isempty(strain_seq)
                # Compute position adjustment via indel DB
                shift = get_indel_shift(ctx.indel_db, annotation.transcript_id, v.strain, annotation.pos_in_cds)
                adjusted_pos = annotation.pos_in_cds + shift

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
end

"""
    build_reference_variation(variations, annotation, seq_id, location, reference_strain)
        -> Variation

Creates reference strain variation record.
Computes adjacent_snp_causes_product_difference.
"""
function build_reference_variation(
    variations::Vector{Variation},
    annotation::PositionAnnotation,
    seq_id::String,
    location::Int,
    reference_strain::String
)
    ref_allele = variations[1].reference
    # If reference field is empty (SNP-only lines), use the base from the first entry
    if isempty(ref_allele)
        ref_allele = variations[1].base
    end

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

    # Mark has_nonsynonomous on each variation
    for v in variations
        if annotation.is_coding == 1
            product_str = join(v.product, ":")
            v.has_nonsynonomous = (product_str != annotation.ref_product) ? 1 : 0
        end
    end

    Variation(
        seq_id,
        location,
        reference_strain,
        ref_allele,
        ref_allele,       # base = reference allele
        "", "",           # coverage, percent (empty for reference)
        "", "",           # quality, pvalue
        "",               # snp_source_id
        annotation.is_coding,
        annotation.pos_in_cds,
        annotation.pos_in_codon_val,
        0,                # downstream_of_frameshift (computed later)
        annotation.transcript_id,
        annotation.is_coding == 1 ? [annotation.ref_product] : String[],
        annotation.ref_codon,
        annotation.ref_codon,        # codon = ref_codon for reference
        adjacent_snp_causes_product_difference,
        annotation.cds_number,
        1                 # matches_reference
    )
end

"""
    fill_coverage_gaps!(variations, annotation, seq_id, location, ctx)

Adds coverage-only variations for missing strains.
Appends to variations vector.
"""
function fill_coverage_gaps!(
    variations::Vector{Variation},
    annotation::PositionAnnotation,
    seq_id::String,
    location::Int,
    ctx::ProcessingContext
)
    # Get reference allele
    ref_allele = variations[1].reference
    if isempty(ref_allele)
        ref_allele = variations[1].base
    end

    # Find adjacent_snp_causes_product_difference from reference variation (last one added)
    adjacent_snp_causes_product_difference = 0
    for v in variations
        if v.strain == ctx.reference_strain
            adjacent_snp_causes_product_difference = v.has_nonsynonomous
            break
        end
    end

    variation_strains = Set(v.strain for v in variations)

    for strain in ctx.all_strains
        strain in variation_strains && continue

        cr = ctx.coverage_readers[strain]
        cov_result = get_coverage(cr, seq_id, location)
        if !isnothing(cov_result)
            cv = Variation()
            cv.sequence_source_id = seq_id
            cv.location = location
            cv.strain = strain
            cv.reference = ref_allele
            cv.base = ref_allele
            cv.coverage = cov_result[1]
            cv.percent = cov_result[2]
            cv.matches_reference = 1
            cv.is_coding = annotation.is_coding
            cv.transcript = annotation.transcript_id
            cv.position_in_cds = annotation.pos_in_cds
            cv.position_in_codon = annotation.pos_in_codon_val
            cv.reference_codon = annotation.ref_codon
            cv.codon = annotation.ref_codon
            cv.product = annotation.is_coding == 1 ? [annotation.ref_product] : String[]
            cv.has_nonsynonomous = adjacent_snp_causes_product_difference
            cv.cds_number = annotation.cds_number
            push!(variations, cv)
        end
    end
end

# ---------------------------------------------------------------------------
# Output writing functions
# ---------------------------------------------------------------------------

"""
    write_cache_entries(cache_fh, variations, ctx)

Computes downstream_of_frameshift and writes cache lines for all non-reference variations.
"""
function write_cache_entries(
    cache_fh::IOStream,
    variations::Vector{Variation},
    ctx::ProcessingContext
)
    # Get reference allele for matches_reference computation
    ref_allele = variations[1].reference
    if isempty(ref_allele)
        ref_allele = variations[1].base
    end

    for v in variations
        if !isempty(v.transcript)
            v.downstream_of_frameshift = check_downstream_of_frameshift(
                ctx.fs_info, v.strain, v.transcript, v.location)
        else
            v.downstream_of_frameshift = 0
        end

        # Skip reference strain in cache output
        v.strain == ctx.reference_strain && continue

        # matches_reference
        v.matches_reference = (v.base == ref_allele) ? 1 : 0

        # Write cache line (20 columns)
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
end

"""
    write_snp_feature(snp_fh, variations, annotation, seq_id, location, reference_strain)

Aggregates statistics and writes SNP feature summary line.
"""
function write_snp_feature(
    snp_fh::IOStream,
    variations::Vector{Variation},
    annotation::PositionAnnotation,
    seq_id::String,
    location::Int,
    reference_strain::String
)
    # Get reference allele
    ref_allele = variations[1].reference
    if isempty(ref_allele)
        ref_allele = variations[1].base
    end

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

"""
    write_allele_and_product_files(allele_fh, product_fh, variations, annotation)

Writes allele statistics and product details for coding variants.
"""
function write_allele_and_product_files(
    allele_fh::IOStream,
    product_fh::IOStream,
    variations::Vector{Variation},
    annotation::PositionAnnotation
)
    annotation.is_coding != 1 && return

    # Allele file: allele, distinct_strain_count, allele_count, average_coverage, average_read_percent
    allele_counts = Dict{String,Int}()
    for v in variations
        allele_counts[v.base] = get(allele_counts, v.base, 0) + 1
    end

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
# Main loop orchestration
# ---------------------------------------------------------------------------

"""
    process_single_position!(cache_pf, snp_pf, seq_id, location, ctx, writers, transcript_cache) -> Bool

Process a single genomic position through the full pipeline.
Returns true if position was processed, false if skipped.
"""
function process_single_position!(
    cache_pf::PeekedFile,
    snp_pf::PeekedFile,
    seq_id::String,
    location::Int,
    ctx::ProcessingContext,
    writers::OutputWriters,
    transcript_cache::TranscriptSequenceCache
)
    # Collect all variations at this position
    (variations, _) = collect_variations_at_position(
        cache_pf, snp_pf, seq_id, location, ctx.undone_strains
    )

    # Skip if empty batch
    isempty(variations) && return false

    # Annotate position
    annotation = annotate_position(seq_id, location, ctx, transcript_cache)

    # Annotate each variation with coding info
    annotate_variations!(variations, annotation, ctx, transcript_cache)

    # Build reference variation record
    ref_variation = build_reference_variation(
        variations, annotation, seq_id, location, ctx.reference_strain
    )

    # Add reference variation to the batch
    push!(variations, ref_variation)

    # Check if there is actual variation (more than one distinct allele)
    if !has_variation(variations)
        return false  # no variation, skip
    end

    # Fill coverage for strains not in the variation batch
    fill_coverage_gaps!(variations, annotation, seq_id, location, ctx)

    # Write outputs
    write_cache_entries(writers.cache_fh, variations, ctx)
    write_snp_feature(writers.snp_fh, variations, annotation, seq_id, location, ctx.reference_strain)
    write_allele_and_product_files(writers.allele_fh, writers.product_fh, variations, annotation)

    true
end

"""
    process_all_positions!(cache_pf, snp_pf, ctx, writers, transcript_cache) -> Int

Main processing loop: sorted merge of cache and SNP files.
Returns count of processed positions.
"""
function process_all_positions!(
    cache_pf::PeekedFile,
    snp_pf::PeekedFile,
    ctx::ProcessingContext,
    writers::OutputWriters,
    transcript_cache::TranscriptSequenceCache
)
    count = 0

    while !cache_pf.exhausted || !snp_pf.exhausted
        # Determine the next (seq_id, location) to process
        next_pos = determine_next_position(cache_pf, snp_pf)
        isnothing(next_pos) && break

        current_seq_id, current_location = next_pos

        # Process this position
        if process_single_position!(
            cache_pf, snp_pf, current_seq_id, current_location,
            ctx, writers, transcript_cache
        )
            count += 1
            if count % 1000 == 0
                @info "Processed $count SNPs"
            end
        end
    end

    count
end

# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------

function main()
    # 1. Parse arguments
    args = parse_args(ARGS)

    # 2. Initialize processing context
    ctx = initialize_processing_context(args)

    # 3. Open output writers
    (writers, temp_cache_file) = open_output_writers(args["cache_file"])

    # 4. Open input files with peek
    cache_pf = open_peeked(args["cache_file"])
    snp_pf = open_peeked(args["snp_file"])

    # 5. Initialize transcript sequence cache
    transcript_cache = TranscriptSequenceCache("", Dict{String,String}())

    # 6. Main processing loop
    process_all_positions!(
        cache_pf, snp_pf, ctx, writers, transcript_cache
    )

    # 7. Close input files
    close_peeked(cache_pf)
    close_peeked(snp_pf)

    # 8. Close output writers
    close_output_writers(writers)

    # 9. Close processing context
    close_processing_context(ctx)

    # 10. Finalize files
    finalize_output_files(
        args["cache_file"], temp_cache_file,
        args["snp_file"], args["undone_strains_file"]
    )
end

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
main()
