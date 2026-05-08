# Plan: Rewrite processSequenceVariations in Julia

## Scope
Rewrite `bin/processSequenceVariationsNew.pl` as `bin/processSequenceVariations.jl`. The upstream preparation of CDS sequences and CDS-space indels (from genomic indels + GTF + consensus FASTA) is now in scope as separate Nextflow processes. Update the Nextflow process definition and Dockerfile accordingly.

## Architecture Overview

### Sequence access — single SQLite DB:
All CDS sequences (reference AND per-strain) live in one SQLite database. No genomic FASTA access needed at all — the upstream process splices CDS sequences for every strain including the reference.

**Coding sequences DB** (`codingSequences.db`) — produced by the upstream CDS-prep process. Schema:
```sql
CREATE TABLE coding_sequences (
  strain TEXT NOT NULL,
  transcript_id TEXT NOT NULL,
  sequence TEXT NOT NULL,
  PRIMARY KEY (strain, transcript_id)
);
```
- The reference strain is included as a row (e.g. `strain='pfal3D7'`).
- When the script moves to a new CDS, it fetches ALL strains for that transcript_id in a single query: `SELECT strain, sequence FROM coding_sequences WHERE transcript_id=?`. Result is cached as a `Dict{String, String}` (strain → sequence). Individual codon lookups then hit this Dict — O(1), no SQL.
- Cache is cleared when moving to the next CDS — keeps memory bounded.
- One query per CDS regardless of strain count. Scales to 1000s of genomes efficiently.

**Indels DB** (`codingIndels.db`) — produced by the upstream CDS-prep process. Schema:
```sql
CREATE TABLE indels (
  strain TEXT NOT NULL,
  transcript_id TEXT NOT NULL,
  position INTEGER NOT NULL,    -- position in CDS, reference coordinates
  shift_amount INTEGER NOT NULL -- per-indel shift (not cumulative; e.g. +3 for insertion, -1 for deletion)
);
CREATE INDEX idx_indels ON indels(transcript_id, strain, position);
```
- **Position adjustment query** (called per strain per coding variant):
  ```sql
  SELECT COALESCE(SUM(shift_amount), 0) FROM indels
  WHERE transcript_id = ? AND strain = ? AND position < ?
  ```
  Returns the total shift to apply to reference CDS position P → adjusted position P' in the per-strain CDS sequence.
- **Frameshift detection** (computed once at startup): for each (strain, CDS) that has indels, fetch all indels sorted by position, compute running cumulative sum, and record whether any prefix sum is not divisible by 3 (and if so, at what position). Stored in a `Dict` in memory for O(1) lookup during the main loop.

**Position adjustment for per-strain sequences:** Per-strain CDS sequences incorporate indels relative to the reference, so a variant at reference CDS position P maps to a different position P' in the per-strain sequence. The rule is simple: P' = P + (sum of all indel shifts in this CDS for this strain where indel_position < P). This sum is computed via a single SQL query against the indels DB (see indels DB schema below). Used only internally for codon extraction — not an output column.

### Why this fixes the performance:
- The Perl `getVariations` uses `sed "${line}q;d"` per line → O(N²) subprocess spawns. **Replaced with sequential file reads** (two open handles, single-pass merge).
- The Perl `getCodingSequence` shells out to `samtools faidx` per call. **Replaced with one bulk SQLite SELECT per CDS** (fetches all strains at once) + in-memory Dict cache. Zero subprocesses. The inner loop is pure Dict lookups + string slicing.
- The Perl `addTranscript` linear-scans all transcripts per variant. **Replaced with sorted interval array + binary search** → O(log N).
- The Perl coordinate-shifting logic (200+ lines) for per-strain CDS extraction is **eliminated** — CDS sequences are pre-spliced and pre-shifted by the upstream process.

---

## Files to Modify

| File | Change |
|---|---|
| `bin/processSequenceVariations.jl` | **NEW** — the Julia script (replaces the Perl script's role) |
| `modules/mergeExperiments.nf` | Add `makeGenomicIndelDb`, `makeCodingData`, `makeCodingIndelsDb`, `makeCodingSequencesDb` processes; update `processSeqVars` inputs |
| `workflows/mergeExperiments.nf` | Wire new CDS-prep processes; update `processSeqVars` call site |
| `Dockerfile` | Add Julia runtime + SQLite.jl precompilation |

### New upstream processes (CDS-prep pipeline)

| Process | Inputs | Outputs | Description |
|---|---|---|---|
| `makeGenomicIndelDb` | `indels.tsv` (combined) | `genomicIndels.db` | Load per-strain genomic indels into SQLite |
| `makeCodingData` | `genomicIndels.db`, GTF, consensus FASTAs | `codingSequences.fa`, `codingIndels.tsv` | Splice CDS sequences per strain; project genomic indels to CDS coordinates |
| `makeCodingSequencesDb` | `codingSequences.fa` | `codingSequences.db` | Load spliced CDS sequences into SQLite |
| `makeCodingIndelsDb` | `codingIndels.tsv` | `codingIndels.db` | Load CDS-space indels into SQLite with index |

The Perl script (`bin/processSequenceVariationsNew.pl`) is kept in place — it's still copied into the image by `ADD /bin/*`. It just won't be called anymore. Remove it in a follow-up cleanup if desired.

---

## Inputs to the Julia Script

| Input | Source | Notes |
|---|---|---|
| `--snp_file` | `snpFile` (existing) | 10-col TSV: seqId, location, strain, ref, base, coverage, percent, quality, pvalue, snp_source_id |
| `--cache_file` | `cacheFile` (existing) | 20-col TSV, same columns as output cache. May be empty on first run. |
| `--undone_strains_file` | `undoneStrainsFile` (existing) | One strain name per line |
| `--coverage_directory` | `coverageDir` (existing) | Per-strain `*.coverage.txt` files |
| `--cds_db` | **NEW** | Path to SQLite DB (`codingSequences.db`) with CDS sequences for all strains including reference |
| `--indel_db` | **NEW** | SQLite DB (`codingIndels.db`) with per-strain indels in CDS coordinates. Used for position adjustment and frameshift detection. See schema below. |
| `--gtf_file` | `gtfFile` (existing) | Gene annotations (CDS intervals for transcript/exon lookup) |
| `--reference_strain` | param (existing) | Strain name to use when looking up reference sequences in transcript DB |

## Outputs (unchanged format)

All four files must match column order exactly. **Note the typo on column 19 of the cache: `has_nonsynonomous` (missing 'y') — must be preserved.**

### Cache / variationFeature.dat — 20 columns
`sequence_source_id, location, strain, reference, base, coverage, percent, quality, pvalue, snp_source_id, is_coding, position_in_cds, position_in_codon, downstream_of_frameshift, transcript, product, reference_codon, codon, has_nonsynonomous, cds_number`

**Removed vs. Perl output:** `shifted_location` and `current_shift` (no longer needed — coordinate shifting is handled upstream).

`product` is colon-joined (e.g. `"L:I"`).

### SNP feature (snpFeature.dat → renamed to variationFeature.dat by Nextflow) — 18 columns
`location, transcript_id, source_id, reference_strain, reference_na, has_nonsynonymous_allele, major_allele, minor_allele, major_allele_count, minor_allele_count, major_product, minor_product, distinct_strain_count, distinct_allele_count, has_coding_mutation, total_allele_count, has_stop_codon, ref_codon`

Major/minor sorted: descending by count, then ascending alphabetically as tiebreaker.

### Allele file — 5 columns
`allele, distinct_strain_count, allele_count, average_coverage, average_read_percent`

Averages formatted to 2 decimal places.

### Product file — 8 columns
`codon, position_in_protein, transcript, count, product, ref_location_cds, ref_location_protein, downstream_of_frameshift`

One row per expanded (non-ambiguous) codon per variant position. IUPAC ambiguity codes in the codon are expanded to all possible codons via Cartesian product.

---

## Julia Script Structure

### 1. Startup / initialization
- Parse CLI args (use a simple hand-rolled parser — no package dependency needed for ~10 flags)
- Open SQLite DB connections: `codingSequences.db` (CDS sequences) and `codingIndels.db` (CDS-space indels)
- Parse GTF → build sorted `CDSInterval` array + per-CDS exon list (exon order and lengths needed for position_in_cds computation)
- Precompute frameshift info from indels DB: query all (strain, transcript_id) combinations that have indels, compute running cumulative sum per combination, record first position where cumulative sum % 3 ≠ 0. Store in `frameshift_info::Dict{String, Dict{String, Tuple{Bool, Int}}}` (strain → transcript_id → (has_frameshift, frameshift_position))
- Open all per-strain coverage files; build `Dict{String, CoverageReader}`
- Open output file handles (cache, snpFeature, allele, product)

### 2. Main loop — sorted merge of cache + SNP file
Both files are sorted by (sequence_source_id, location). Single-pass merge using two open file handles with "peeked" next lines. For each genomic position:

```
read all cache rows at this position (skip strains in undone set)
read all SNP rows at this position (skip strains already seen from cache)
→ yields one batch of Variation records per position
```

This replaces the O(N²) `sed`-based line reading. Total I/O: one pass over each file.

### 3. Per-position processing pipeline
For each batch of variations at a genomic position:

```
a) Determine (seq_id, location) from the batch
b) Binary search CDSInterval array → is it coding? Which transcript? Which CDS exon?
c) If coding:
     - Compute position_in_cds (sum prior exon lengths + offset in current exon, from GTF)
     - If CDS changed: bulk-fetch all strains' sequences for this transcript_id (one SQL query → Dict cache)
     - Extract reference codon at position_in_cds from cached reference sequence; translate → reference product
     - For each strain variation:
         - Look up strain's sequence from the Dict cache (O(1))
         - Compute adjusted position: position_in_cds + SUM of indels before position_in_cds (SQL query on indels DB)
         - Extract codon at adjusted position
         - Expand IUPAC ambiguities → all possible codons
         - Translate each → products list
d) Compute downstream_of_frameshift from frameshift_info (only for coding variants)
e) Build reference variation record, add to batch
f) Fill in coverage for strains not in the batch (sequential scan of coverage files)
g) Skip position if no actual variation (all strains match reference) — same as Perl
h) Write cache records for all non-reference strains
i) Write SNP feature summary record
j) If coding: write allele and product summary records
```

### 4. Key algorithms

**Sequence access via SQLite + in-memory caching:**
```julia
# When CDS (transcript_id) changes, bulk-fetch all strains at once:
function load_cds_sequences(db, transcript_id) -> Dict{String, String}
    # SELECT strain, sequence FROM coding_sequences WHERE transcript_id = ?
    # Returns Dict mapping strain -> sequence for all strains
end

# Cache: current_cds_seqs::Dict{String, String}
# Cleared and reloaded when transcript_id changes.
# Per-variant lookup is just: current_cds_seqs[strain] → O(1)

# Position adjustment for per-strain sequences:
# For a variant at reference CDS position P in transcript_id T for strain S:
#   shift = SQL: SELECT COALESCE(SUM(shift_amount), 0) FROM indels
#                WHERE transcript_id=T AND strain=S AND position < P
#   adjusted_P = P + shift
#
# adjusted_P indexes into the per-strain CDS sequence.
# For the reference strain there are no indels, so shift=0 and adjusted_P == P.
```

**CDS lookup — binary search on sorted intervals:**
```julia
# intervals sorted by (seq_id, start)
# Binary search for rightmost interval with start <= location on matching seq_id
# Then check location <= end_pos
```
O(log N) per variant. N = total CDS exons across all transcripts.

**Frameshift detection:**
At startup, for each (strain, transcript_id): scan CDS-space indels in order, accumulate shift. Record the position of the first indel where cumulative_shift % 3 ≠ 0. At runtime: if variant's position_in_cds > that position → downstream_of_frameshift = 1. (CDS coordinates are always 5'→3', so comparison direction is always `>` regardless of strand — simplification over the Perl code.)

**IUPAC codon expansion:**
Cartesian product of expanded bases at each position. E.g. codon "RYA" → {A,G} × {C,T} × {A} = ["ACA","ATA","GCA","GTA"]. Each expanded codon is translated independently.

**Coverage file reading:**
Same sequential-scan approach as the Perl code. Each coverage file is sorted by (seq_id, start). Maintain a cursor per strain. For each variant position, advance each strain's cursor until the current span covers the position (or passes it). Extract coverage/percent from the span's arrays by offset.

### 5. Reverse complement
All CDS sequences in the DB are assumed to already be in sense (5'→3') orientation — the upstream CDS-prep process is responsible for reverse-complementing minus-strand genes when splicing CDS. The Julia script does not need to reverse-complement anything. (If this assumption is violated, the upstream process has a bug.)

---

## Nextflow Changes

### modules/mergeExperiments.nf — `processSeqVars` process

**New inputs to add:**
- `path cdsDb` — SQLite DB (`codingSequences.db`) with CDS sequences (all strains + reference)
- `path indelDb` — SQLite DB (`codingIndels.db`) with per-strain indels in CDS coordinates (position adjustment + frameshift detection)

**Remove:** `path genomeFasta`, `path consensusFasta`, and the `gunzip` step (all sequence access now via SQLite)

**Script block:** call `julia /usr/bin/processSequenceVariations.jl` with the new flag set (see Inputs table above). Keep the `mv snpFeature.dat variationFeature.dat` at the end.

### workflows/mergeExperiments.nf — call site

Wire `cdsDb` and `indelDb` to the outputs of the upstream CDS-prep processes (`makeCodingSequencesDb` and `makeCodingIndelsDb`). Update the `processSeqVars(...)` argument list accordingly.

---

## Dockerfile Changes

After the main `apt-get` block, add:
```dockerfile
ENV JULIA_VERSION=1.10.8
RUN wget -q https://julialang-releases.github.io/pub/julia/${JULIA_VERSION}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz \
    && tar xzf julia-${JULIA_VERSION}-linux-x86_64.tar.gz \
    && mv julia-${JULIA_VERSION} /opt/julia \
    && rm julia-${JULIA_VERSION}-linux-x86_64.tar.gz
ENV PATH=/opt/julia/bin:$PATH

# Precompile SQLite.jl (only external dependency)
RUN julia -e 'using Pkg; Pkg.add("SQLite"); using SQLite'
```

The existing `ADD /bin/* /usr/bin/` picks up the new `.jl` file automatically.

---

## Verification

1. **Unit-level:** Write a small test script that exercises each key function (SQLite sequence fetch + position adjustment, CDS binary search, frameshift detection, IUPAC expansion, codon translation, coverage cursor advance) against known inputs/outputs.
2. **Integration:** Run with the existing test data (if any in the repo). Compare all four output files byte-for-byte against the Perl script's outputs on the same inputs. Pay special attention to column order and the `has_nonsynonomous` typo.
3. **Pipeline:** Run `nextflow run main.nf -profile mergeExperiments` end-to-end with the updated process. Confirm outputs land in the expected directories.
