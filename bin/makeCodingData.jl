#!/usr/bin/env julia

# makeCodingData.jl
#
# Inputs (staged flat into work dir by Nextflow):
#   *_consensus.fa.gz     — per-strain consensus FASTAs; strain = basename before _consensus.fa.gz
#   --genomic_indel_db    — SQLite DB of genomic-space indels (genomicIndels.db)
#   --gtf_file            — GTF annotation
#   --genome_fasta        — Reference genome FASTA
#   --cds_db_out          — output path for codingSequences.db
#   --indels_db_out       — output path for codingIndels.db
#
# Outputs:
#   codingSequences.db    — SQLite: coding_sequences(strain, transcript_id, sequence)
#   codingIndels.db       — SQLite: indels(strain, transcript_id, position, shift_amount)

include("$(@__DIR__)/GtfUtils.jl")
include("$(@__DIR__)/Benchmark.jl")

using SQLite
using SQLite.DBInterface: execute

# ---------------------------------------------------------------------------
# CLI argument parsing
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
# FASTA reading
# ---------------------------------------------------------------------------

"""
    read_fasta(path) -> Dict{String, String}

Read a (optionally gzipped) FASTA file. Returns seq_id => sequence (uppercase).
"""
function read_fasta(path::String)
    seqs = Dict{String,String}()
    cmd = endswith(path, ".gz") ? `gunzip -c $path` : `cat $path`
    current_id = ""
    buf = IOBuffer()
    for line in eachline(open(cmd))
        if startswith(line, '>')
            if !isempty(current_id)
                seqs[current_id] = String(take!(buf))
            end
            current_id = String(split(line[2:end])[1])
        else
            write(buf, strip(line))
        end
    end
    if !isempty(current_id)
        seqs[current_id] = String(take!(buf))
    end
    seqs
end

"""
    strip_strain_prefix(seqs, strain) -> Dict{String, String}

Consensus FASTA deflines are written as "<strain>_<sequenceId>" so they are
globally unique across strains (e.g. for DB loading / display). This strips the
exact "<strain>_" prefix so the returned dict is keyed by the bare sequenceId,
matching the sequence_ids used in the GTF for CDS coordinate lookups.

Dropping the strain from the *key* is safe only because the strain identity is
retained separately by the caller: `main()` derives it from the filename, passes
it to `process_strain`, and it becomes the `strain` column (part of the primary
key) in both output SQLite tables. Keys lacking the prefix are left unchanged.
"""
function strip_strain_prefix(seqs::Dict{String,String}, strain::String)
    prefix = strain * "_"
    Dict(
        (startswith(k, prefix) ? k[nextind(k, lastindex(prefix)):end] : k) => v
        for (k, v) in seqs
    )
end

# ---------------------------------------------------------------------------
# CDS sequence extraction
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# In-memory genomic-indel index (one strain at a time)
#
# makeCodingData processes a single strain across all transcripts, so instead of
# issuing a SQL query per exon/transcript (fasta_offset + project_indels_to_cds
# fired ~O(strains × transcripts × exons) queries — the makeCodingData
# bottleneck), we pull a strain's genomic indels ONCE and answer both lookups
# from sorted per-sequence arrays via binary search.
# ---------------------------------------------------------------------------

struct GenomicIndelIndex
    # sequence_id -> (positions, shifts, cumshifts) — all ascending by position;
    # cumshifts is the running prefix sum of shifts.
    byseq::Dict{String, Tuple{Vector{Int}, Vector{Int}, Vector{Int}}}
end

"""
    load_strain_genomic_indels(indel_db, strain) -> GenomicIndelIndex

Load all genomic indels for `strain` in one query and build per-sequence sorted
position/shift arrays plus a prefix sum of shifts.
"""
function load_strain_genomic_indels(indel_db::SQLite.DB, strain::String)::GenomicIndelIndex
    byseq = Dict{String, Tuple{Vector{Int}, Vector{Int}, Vector{Int}}}()
    rows = @bench "sql_load_strain_indels" [(String(r[1]), Int(r[2]), Int(r[3])) for r in execute(indel_db,
        "SELECT sequence_id, position, shift FROM genomic_indels
         WHERE strain = ? ORDER BY sequence_id, position",
        [strain])]
    for (seq_id, pos, shift) in rows
        entry = get(byseq, seq_id, nothing)
        if entry === nothing
            entry = (Int[], Int[], Int[])
            byseq[seq_id] = entry
        end
        push!(entry[1], pos)
        push!(entry[2], shift)
    end
    for (_, (positions, shifts, cums)) in byseq
        c = 0
        for s in shifts
            c += s
            push!(cums, c)
        end
    end
    GenomicIndelIndex(byseq)
end

"""
    fasta_offset(idx, seq_id, before_pos) -> Int

Cumulative indel shift on `seq_id` at all positions strictly less than
`before_pos`. Converts a reference genomic coordinate to the corresponding
position in a consensus FASTA that embeds genomic indels. Matches the old SQL
`SUM(shift) WHERE ... AND position < ?` (0 when none) via binary search.
"""
function fasta_offset(idx::GenomicIndelIndex, seq_id::String, before_pos::Int)::Int
    entry = get(idx.byseq, seq_id, nothing)
    entry === nothing && return 0
    positions, _, cums = entry
    i = searchsortedlast(positions, before_pos - 1)  # position < before_pos
    i == 0 ? 0 : cums[i]
end

"""
    extract_cds_sequence(genome_seqs, exon_list, indel_idx) -> String

Splice and return the CDS sequence for a transcript from genome sequences.
`exon_list` must be in 5'→3' order (as returned by parse_gtf).
Reverse-complements minus-strand transcripts using IUPAC-aware complement.
Returns "" if any exon's seq_id is missing from genome_seqs.

For non-reference strains, pass the strain's `indel_idx` (a GenomicIndelIndex) so
exon slice coordinates are adjusted for upstream genomic indels embedded in the
consensus FASTA. Pass `nothing` for the reference (coordinates used as-is).
"""
function extract_cds_sequence(genome_seqs::Dict{String,String}, exon_list::Vector{CdsExon},
                               indel_idx::Union{GenomicIndelIndex,Nothing}=nothing)
    isempty(exon_list) && return ""
    parts = String[]
    for e in exon_list
        seq = get(genome_seqs, e.seq_id, nothing)
        seq === nothing && return ""
        if indel_idx !== nothing
            fstart = e.start + fasta_offset(indel_idx, e.seq_id, e.start)
            fstop  = e.stop  + fasta_offset(indel_idx, e.seq_id, e.stop + 1)
        else
            fstart, fstop = e.start, e.stop
        end
        push!(parts, seq[fstart:fstop])
    end
    if exon_list[1].strand != '-'
        cds = join(parts)
    else
        cds = join(reverse_complement(p) for p in parts)
    end
    cds
end

# ---------------------------------------------------------------------------
# Genomic → CDS indel projection
# ---------------------------------------------------------------------------

"""
    project_indels_to_cds(indel_idx, exon_list) -> Vector{Tuple{Int,Int}}

Find the strain's genomic indels (from `indel_idx`) that fall within the exons of
a transcript, convert each to a 0-based CDS-space position, and return
[(cds_position, shift_amount), ...] sorted in 5'→3' order.
"""
function project_indels_to_cds(indel_idx::GenomicIndelIndex, exon_list::Vector{CdsExon})
    result = Tuple{Int,Int}[]
    isempty(exon_list) && return result

    seq_id     = exon_list[1].seq_id
    cds_offset = 0

    entry = get(indel_idx.byseq, seq_id, nothing)
    positions, shifts = entry === nothing ? (Int[], Int[]) : (entry[1], entry[2])

    for e in exon_list
        # In-memory equivalent of the old
        #   WHERE position >= e.start AND position <= e.stop ORDER BY position
        lo = searchsortedfirst(positions, e.start)
        hi = searchsortedlast(positions, e.stop)
        for i in lo:hi
            gpos  = positions[i]
            shift = shifts[i]
            cds_pos = if e.strand == '-'
                cds_offset + (e.stop - gpos)
            else
                cds_offset + (gpos - e.start)
            end
            push!(result, (cds_pos, shift))
        end
        cds_offset += e.stop - e.start + 1
    end

    # Minus-strand exons were iterated 3'→5' (descending genomic); reverse to 5'→3'
    if exon_list[1].strand == '-'
        reverse!(result)
    end

    result
end

# ---------------------------------------------------------------------------
# SQLite output helpers
# ---------------------------------------------------------------------------

function create_cds_db(path::String)
    db = SQLite.DB(path)
    execute(db, """
        CREATE TABLE coding_sequences (
          strain        TEXT NOT NULL,
          transcript_id TEXT NOT NULL,
          sequence      TEXT NOT NULL,
          PRIMARY KEY (strain, transcript_id)
        )
    """)
    db
end

function create_indels_db(path::String)
    db = SQLite.DB(path)
    execute(db, """
        CREATE TABLE indels (
          strain        TEXT    NOT NULL,
          transcript_id TEXT    NOT NULL,
          position      INTEGER NOT NULL,
          shift_amount  INTEGER NOT NULL
        )
    """)
    db
end

function finalize_indels_db(db::SQLite.DB)
    @bench "finalize_indels_index" execute(db, "CREATE INDEX idx_indels ON indels(transcript_id, strain, position)")
end

function finalize_cds_db(db::SQLite.DB)
    # processSequenceVariations.jl queries coding_sequences by transcript_id alone.
    # The PRIMARY KEY (strain, transcript_id) autoindex leads with strain, so a
    # transcript_id-only lookup cannot use it and full-scans the (large) table.
    # This index makes that lookup a seek.
    @bench "finalize_cds_index" execute(db, "CREATE INDEX idx_coding_sequences_transcript ON coding_sequences(transcript_id)")
end

# ---------------------------------------------------------------------------
# Per-strain processing
# ---------------------------------------------------------------------------

function process_strain(strain, genome_seqs, by_transcript, indel_src_db,
                        cds_insert_stmt, indels_insert_stmt)
    # One query per strain instead of per (transcript, exon). The reference gets no
    # coordinate offset in extract_cds (indel_idx=nothing), matching prior behavior;
    # project_indels still uses the loaded index (empty for the reference).
    indel_idx  = load_strain_genomic_indels(indel_src_db, strain)
    extract_idx = strain == "reference" ? nothing : indel_idx
    for (transcript_id, exon_list) in by_transcript
        seq = @bench "extract_cds" extract_cds_sequence(genome_seqs, exon_list, extract_idx)
        isempty(seq) && continue
        @bench "cds_insert" execute(cds_insert_stmt, [strain, transcript_id, seq])

        cds_indels = @bench "project_indels" project_indels_to_cds(indel_idx, exon_list)
        for (pos, shift) in cds_indels
            @bench "indels_insert" execute(indels_insert_stmt, [strain, transcript_id, pos, shift])
        end
    end
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    args = parse_args(ARGS)
    global BENCHMARK = haskey(args, "benchmark")
    bench_t0 = time_ns()

    genomic_indel_db = args["genomic_indel_db"]
    gtf_file         = args["gtf_file"]
    genome_fasta     = args["genome_fasta"]
    cds_db_out       = args["cds_db_out"]
    indels_db_out    = args["indels_db_out"]

    println(stderr, "Parsing GTF: $gtf_file")
    _, by_transcript = @bench "parse_gtf" parse_gtf(gtf_file)

    println(stderr, "Reading reference genome: $genome_fasta")
    ref_seqs = @bench "read_fasta" read_fasta(genome_fasta)

    println(stderr, "Opening genomic indels DB: $genomic_indel_db")
    indel_src_db = SQLite.DB(genomic_indel_db)

    cds_db    = create_cds_db(cds_db_out)
    indels_db = create_indels_db(indels_db_out)

    cds_insert_stmt    = SQLite.Stmt(cds_db,
        "INSERT INTO coding_sequences(strain, transcript_id, sequence) VALUES (?, ?, ?)")
    indels_insert_stmt = SQLite.Stmt(indels_db,
        "INSERT INTO indels(strain, transcript_id, position, shift_amount) VALUES (?, ?, ?, ?)")

    # Reference strain — sequences from genomeFastaFile, no indels
    println(stderr, "Extracting reference CDS sequences")
    SQLite.transaction(cds_db) do
        process_strain("reference", ref_seqs, by_transcript, indel_src_db,
                       cds_insert_stmt, indels_insert_stmt)
    end

    # Per-strain consensus FASTAs staged flat in work dir
    fasta_files = sort(filter(f -> endswith(f, "_consensus.fa.gz"), readdir(".")))
    for fasta_path in fasta_files
        strain = replace(basename(fasta_path), "_consensus.fa.gz" => "")
        println(stderr, "Processing strain: $strain")
        strain_seqs = strip_strain_prefix(@bench("read_fasta", read_fasta(fasta_path)), strain)
        SQLite.transaction(cds_db) do
            SQLite.transaction(indels_db) do
                process_strain(strain, strain_seqs, by_transcript, indel_src_db,
                               cds_insert_stmt, indels_insert_stmt)
            end
        end
    end

    finalize_indels_db(indels_db)
    finalize_cds_db(cds_db)
    println(stderr, "Done. Wrote $cds_db_out and $indels_db_out")

    if BENCHMARK
        n_strains = 1 + length(fasta_files)   # reference + per-strain consensus FASTAs
        print_benchmark_report(stderr, time_ns() - bench_t0;
            summary = "strains (incl. reference): $n_strains   transcripts: $(length(by_transcript))",
            per_label = "calls/strain", per_denom = n_strains)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
