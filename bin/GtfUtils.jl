#!/usr/bin/env julia

# GtfUtils.jl
# Shared GTF parsing, CDS interval, and sequence utilities.
# Included by makeCodingData.jl and processSequenceVariations.jl via:
#   include("$(@__DIR__)/GtfUtils.jl")

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

struct CdsExon
    seq_id        :: String
    start         :: Int      # 1-based, inclusive (genomic)
    stop          :: Int      # 1-based, inclusive (genomic)
    strand        :: Char     # '+' or '-'
    transcript_id :: String
    exon_number   :: Int      # 1-based order in transcript (5'→3')
end

# Sorted array of all CDS exons across all transcripts, ordered by
# (seq_id, start) for binary search during variant lookup.
const CDSInterval = CdsExon

# Per-transcript ordered list of exons (5'→3')
# transcript_id => Vector{CdsExon} sorted in transcription order
const TranscriptExons = Dict{String, Vector{CdsExon}}

# ---------------------------------------------------------------------------
# IUPAC reverse complement
# ---------------------------------------------------------------------------

const IUPAC_COMPLEMENT = Dict(
    'A'=>'T', 'T'=>'A', 'G'=>'C', 'C'=>'G',
    'R'=>'Y', 'Y'=>'R', 'S'=>'S', 'W'=>'W',
    'K'=>'M', 'M'=>'K', 'B'=>'V', 'V'=>'B',
    'D'=>'H', 'H'=>'D', 'N'=>'N',
    'a'=>'t', 't'=>'a', 'g'=>'c', 'c'=>'g',
    'r'=>'y', 'y'=>'r', 's'=>'s', 'w'=>'w',
    'k'=>'m', 'm'=>'k', 'b'=>'v', 'v'=>'b',
    'd'=>'h', 'h'=>'d', 'n'=>'n',
)

"""
    reverse_complement(seq) -> String

Return the reverse complement of `seq` using full IUPAC ambiguity codes.
"""
function reverse_complement(seq::String)
    String([get(IUPAC_COMPLEMENT, c, 'N') for c in reverse(seq)])
end

# ---------------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------------

"""
    parse_gtf(gtf_file) -> (Vector{CDSInterval}, TranscriptExons)

Parse a GTF file and return:
- `intervals`: all CDS exons sorted by (seq_id, start) for binary search
- `by_transcript`: per-transcript exon lists sorted in 5'→3' order
"""
function parse_gtf(gtf_file::String)
    exons = CdsExon[]

    open(gtf_file) do fh
        for line in eachline(fh)
            startswith(line, '#') && continue
            fields = split(line, '\t')
            length(fields) < 9 && continue
            fields[3] == "CDS" || continue

            seq_id = fields[1]
            start  = parse(Int, fields[4])
            stop   = parse(Int, fields[5])
            strand = fields[7][1]
            attrs  = fields[9]

            transcript_id = _extract_attr(attrs, "transcript_id")
            isempty(transcript_id) && continue

            exon_number_str = _extract_attr(attrs, "exon_number")
            exon_number = isempty(exon_number_str) ? 0 : parse(Int, exon_number_str)

            push!(exons, CdsExon(seq_id, start, stop, strand, transcript_id, exon_number))
        end
    end

    # Sort intervals by (seq_id, start) for binary search
    intervals = sort(exons, by = e -> (e.seq_id, e.start))

    # Group by transcript and sort in 5'→3' order
    by_transcript = TranscriptExons()
    for e in exons
        push!(get!(by_transcript, e.transcript_id, CdsExon[]), e)
    end
    for (_, exon_list) in by_transcript
        if !isempty(exon_list) && exon_list[1].strand == '-'
            sort!(exon_list, by = e -> e.start, rev = true)
        else
            sort!(exon_list, by = e -> e.start)
        end
        # Re-assign exon_number based on 5'→3' position order if absent from GTF
        if all(e -> e.exon_number == 0, exon_list)
            for i in eachindex(exon_list)
                e = exon_list[i]
                exon_list[i] = CdsExon(e.seq_id, e.start, e.stop, e.strand, e.transcript_id, i)
            end
        end
    end

    return intervals, by_transcript
end

function _extract_attr(attrs::AbstractString, key::String)
    pattern = Regex("(?:^|;)\\s*" * key * "\\s+\"([^\"]+)\"")
    m = match(pattern, attrs)
    m === nothing ? "" : String(m.captures[1])
end

# ---------------------------------------------------------------------------
# CDS lookup — binary search on sorted intervals
# ---------------------------------------------------------------------------

"""
    find_cds_exon(intervals, seq_id, location) -> Union{CdsExon, Nothing}

Binary search the sorted interval array for the CDS exon covering `location`
on `seq_id`. Returns `nothing` if not in any CDS.
"""
function find_cds_exon(intervals::Vector{CDSInterval}, seq_id::String, location::Int)
    lo, hi = 1, length(intervals)
    result = nothing
    while lo <= hi
        mid = (lo + hi) ÷ 2
        e = intervals[mid]
        if e.seq_id < seq_id || (e.seq_id == seq_id && e.start <= location)
            if e.seq_id == seq_id && e.start <= location <= e.stop
                result = e
            end
            lo = mid + 1
        else
            hi = mid - 1
        end
    end
    return result
end

# ---------------------------------------------------------------------------
# CDS position calculation
# ---------------------------------------------------------------------------

"""
    position_in_cds(exon, location, by_transcript) -> Int

Compute the 0-based position of `location` (genomic, 1-based) within the
spliced CDS sequence for `exon`'s transcript (5'→3' orientation).
"""
function position_in_cds(exon::CdsExon, location::Int, by_transcript::TranscriptExons)
    exon_list = by_transcript[exon.transcript_id]
    pos = 0
    for e in exon_list
        e.exon_number == exon.exon_number && break
        pos += e.stop - e.start + 1
    end
    if exon.strand == '-'
        pos += exon.stop - location
    else
        pos += location - exon.start
    end
    return pos
end
