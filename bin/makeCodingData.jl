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

"""
    fasta_offset(indel_db, strain, seq_id, before_pos) -> Int

Return the cumulative indel shift for `strain` on `seq_id` at all positions
strictly less than `before_pos`. Used to convert a reference genomic coordinate
to the corresponding position in a consensus FASTA that embeds genomic indels.
"""
function fasta_offset(indel_db::SQLite.DB, strain::String, seq_id::String, before_pos::Int)
    row = first(execute(indel_db,
        "SELECT COALESCE(SUM(shift), 0) FROM genomic_indels
         WHERE strain = ? AND sequence_id = ? AND position < ?",
        [strain, seq_id, before_pos]))
    row[1]
end

"""
    extract_cds_sequence(genome_seqs, exon_list, strain, indel_db) -> String

Splice and return the CDS sequence for a transcript from genome sequences.
`exon_list` must be in 5'→3' order (as returned by parse_gtf).
Reverse-complements minus-strand transcripts using IUPAC-aware complement.
Returns "" if any exon's seq_id is missing from genome_seqs.

For non-reference strains, pass `strain` and `indel_db` so that exon slice
coordinates are adjusted for upstream genomic indels embedded in the consensus
FASTA.
"""
function extract_cds_sequence(genome_seqs::Dict{String,String}, exon_list::Vector{CdsExon},
                               strain::String="reference",
                               indel_db::Union{SQLite.DB,Nothing}=nothing)
    isempty(exon_list) && return ""
    parts = String[]
    for e in exon_list
        seq = get(genome_seqs, e.seq_id, nothing)
        seq === nothing && return ""
        if indel_db !== nothing
            fstart = e.start + fasta_offset(indel_db, strain, e.seq_id, e.start)
            fstop  = e.stop  + fasta_offset(indel_db, strain, e.seq_id, e.stop + 1)
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
    project_indels_to_cds(db, strain, exon_list) -> Vector{Tuple{Int,Int}}

Query genomic indels for `strain` that fall within the exons of a transcript,
convert each to a 0-based CDS-space position, and return
[(cds_position, shift_amount), ...] sorted in 5'→3' order.
"""
function project_indels_to_cds(db::SQLite.DB, strain::String, exon_list::Vector{CdsExon})
    result = Tuple{Int,Int}[]
    isempty(exon_list) && return result

    seq_id     = exon_list[1].seq_id
    cds_offset = 0

    for e in exon_list
        rows = execute(db,
            "SELECT position, shift FROM genomic_indels
             WHERE strain = ? AND sequence_id = ? AND position >= ? AND position <= ?
             ORDER BY position",
            [strain, seq_id, e.start, e.stop])
        for row in rows
            gpos  = row[1]
            shift = row[2]
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
    execute(db, "CREATE INDEX idx_indels ON indels(transcript_id, strain, position)")
end

function finalize_cds_db(db::SQLite.DB)
    # processSequenceVariations.jl queries coding_sequences by transcript_id alone.
    # The PRIMARY KEY (strain, transcript_id) autoindex leads with strain, so a
    # transcript_id-only lookup cannot use it and full-scans the (large) table.
    # This index makes that lookup a seek.
    execute(db, "CREATE INDEX idx_coding_sequences_transcript ON coding_sequences(transcript_id)")
end

# ---------------------------------------------------------------------------
# Per-strain processing
# ---------------------------------------------------------------------------

function process_strain(strain, genome_seqs, by_transcript, indel_src_db,
                        cds_insert_stmt, indels_insert_stmt)
    db_for_coords = strain == "reference" ? nothing : indel_src_db
    for (transcript_id, exon_list) in by_transcript
        seq = extract_cds_sequence(genome_seqs, exon_list, strain, db_for_coords)
        isempty(seq) && continue
        execute(cds_insert_stmt, [strain, transcript_id, seq])

        cds_indels = project_indels_to_cds(indel_src_db, strain, exon_list)
        for (pos, shift) in cds_indels
            execute(indels_insert_stmt, [strain, transcript_id, pos, shift])
        end
    end
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    args = parse_args(ARGS)

    genomic_indel_db = args["genomic_indel_db"]
    gtf_file         = args["gtf_file"]
    genome_fasta     = args["genome_fasta"]
    cds_db_out       = args["cds_db_out"]
    indels_db_out    = args["indels_db_out"]

    println(stderr, "Parsing GTF: $gtf_file")
    _, by_transcript = parse_gtf(gtf_file)

    println(stderr, "Reading reference genome: $genome_fasta")
    ref_seqs = read_fasta(genome_fasta)

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
        strain_seqs = strip_strain_prefix(read_fasta(fasta_path), strain)
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
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
