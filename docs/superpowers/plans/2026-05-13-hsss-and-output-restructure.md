# HSSS + Output File Restructure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Strip variationFeature.dat to 14 genomic columns, merge product.dat+cache.tsv into transcript_product.dat (12 cols, per-transcript), and add HSSS binary strain files for 4 read-frequency cutoffs.

**Architecture:** All changes are inside `bin/processSequenceVariations.jl`. `write_snp_feature` moves outside the annotation loop and loses CDS columns. `write_product_file` is renamed to `write_transcript_product`, absorbs the cache write, and gains pos_in_cds/pos_in_protein/downstream_of_frameshift columns. A new `HsssState` struct manages 4×(N_strains+2) file handles written in a single pass. `OutputWriters` loses `cache_fh` and gains `transcript_product_fh` + `hsss`.

**Tech Stack:** Julia 1.10, Nextflow DSL2, Python pytest (e2e), Perl Test2::V0 (not touched here)

---

## File Map

| File | Changes |
|---|---|
| `bin/processSequenceVariations.jl` | All logic changes — 4 tasks |
| `modules/mergeExperiments.nf` | Output wiring — Task 4 |
| `testing/t/handleVariantRecord.jl` | Unit tests — Tasks 1–3 |
| `testing/t/test_mergeExperiments_e2e.py` | E2e expectations — Task 5 |

---

## Task 1: Strip variationFeature.dat to 14 genomic columns

**Files:**
- Modify: `testing/t/handleVariantRecord.jl` (add failing test)
- Modify: `bin/processSequenceVariations.jl:1264-1349` (`write_snp_feature`), `:1078` (header), `:1847-1868` (call site in `handle_variant_record!`)

### New column set (14, 0-indexed for Python tests)

| idx | name |
|---|---|
| 0 | location |
| 1 | seq_id |
| 2 | reference_strain |
| 3 | ref_allele |
| 4 | major_allele |
| 5 | minor_allele |
| 6 | major_allele_strain_count |
| 7 | minor_allele_strain_count |
| 8 | major_allele_frequency |
| 9 | minor_allele_frequency |
| 10 | distinct_strain_count |
| 11 | distinct_allele_count |
| 12 | total_ploidy_count |
| 13 | is_coding |

- [ ] **Step 1: Write failing unit test**

Append to `testing/t/handleVariantRecord.jl`:

```julia
# ---------------------------------------------------------------------------
# write_snp_feature — 14 genomic columns, called once per position
# ---------------------------------------------------------------------------

@testset "write_snp_feature emits 14 columns, no CDS fields" begin
    sample_id_map = Dict{String,Int}("ref" => 1, "s1" => 2, "s2" => 3)

    # ref strain: A; s1: T; s2: A (matches ref)
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "T"; v_s1.reference  = "A"; v_s1.ploidy = 1
    v_s2  = Variation(); v_s2.strain  = "s2";  v_s2.base  = "A"; v_s2.reference  = "A"; v_s2.ploidy = 1

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1, v_s2], 1, "LmjF.01", 500, "ref", sample_id_map)
    lines = filter(!isempty, split(String(take!(buf)), "\n"))

    @test length(lines) == 1
    fields = split(lines[1], "\t")
    @test length(fields) == 14
    @test fields[1]  == "500"       # location
    @test fields[2]  == "LmjF.01"  # seq_id
    @test fields[3]  == "ref"       # reference_strain
    @test fields[4]  == "A"         # ref_allele
    @test fields[14] == "1"         # is_coding
end

@testset "write_snp_feature is_coding=0 for non-coding position" begin
    sample_id_map = Dict{String,Int}("ref" => 1, "s1" => 2)
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"; v_ref.ploidy = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "T"; v_s1.reference  = "A"; v_s1.ploidy = 1

    buf = IOBuffer()
    write_snp_feature(buf, [v_ref, v_s1], 0, "LmjF.01", 200, "ref", sample_id_map)
    fields = split(filter(!isempty, split(String(take!(buf)), "\n"))[1], "\t")
    @test fields[14] == "0"
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia testing/t/handleVariantRecord.jl 2>&1 | tail -20
```

Expected: `MethodError` — `write_snp_feature` signature mismatch (currently takes `annotation`, not `is_coding::Int`).

- [ ] **Step 3: Replace `write_snp_feature` body and signature**

Replace the entire `write_snp_feature` function (lines ~1264–1349) with:

```julia
function write_snp_feature(
    snp_fh::IO,
    variations::Vector{Variation},
    is_coding::Int,
    seq_id::String,
    location::Int,
    reference_strain::String,
    sample_id_map::Dict{String,Int}
)
    ref_allele = variations[1].reference
    isempty(ref_allele) && (ref_allele = variations[1].base)

    allele_counts        = Dict{String,Int}()
    allele_ploidy_counts = Dict{String,Int}()
    strain_set           = Set{String}()
    total_ploidy_count   = 0

    for v in variations
        if !isempty(v.alt_allele)
            allele_counts[v.reference]  = get(allele_counts, v.reference,  0) + 1
            allele_counts[v.alt_allele] = get(allele_counts, v.alt_allele, 0) + 1
            allele_ploidy_counts[v.reference]  = get(allele_ploidy_counts, v.reference,  0) + 1
            allele_ploidy_counts[v.alt_allele] = get(allele_ploidy_counts, v.alt_allele, 0) + 1
            total_ploidy_count += 2
        else
            allele_counts[v.base]        = get(allele_counts, v.base, 0) + 1
            allele_ploidy_counts[v.base] = get(allele_ploidy_counts, v.base, 0) + v.ploidy
            total_ploidy_count += v.ploidy
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
    major_allele_frequency    = @sprintf("%.4f", allele_ploidy_counts[major_allele] / total_ploidy_count)
    minor_allele_frequency    = length(sorted_alleles) > 1 ?
        @sprintf("%.4f", allele_ploidy_counts[minor_allele] / total_ploidy_count) : ""

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
        string(is_coding)
    ], "\t"), "\n")
end
```

- [ ] **Step 4: Update the header line in `open_output_writers` (~line 1078)**

```julia
write(snp_fh, "location\tseq_id\treference_strain\tref_allele\tmajor_allele\tminor_allele\tmajor_allele_strain_count\tminor_allele_strain_count\tmajor_allele_frequency\tminor_allele_frequency\tdistinct_strain_count\tdistinct_allele_count\ttotal_ploidy_count\tis_coding\n")
```

- [ ] **Step 5: Move `write_snp_feature` call outside annotation loop in `handle_variant_record!`**

Before the annotation loop, add:
```julia
first_all_vars = nothing
```

Inside the annotation loop, right after `any_output = true`, add:
```julia
if isnothing(first_all_vars)
    first_all_vars = all_vars
end
```

Remove the existing `write_snp_feature(...)` call from inside the loop.

After the annotation loop, replace `any_output || return false` with:
```julia
any_output || return false

is_coding = any(a.is_coding == 1 for a in annotations) ? 1 : 0
write_snp_feature(writers.snp_fh, first_all_vars, is_coding, seq_id, location,
                  ctx.reference_strain, ctx.sample_id_map)
```

- [ ] **Step 6: Run tests to verify they pass**

```bash
julia testing/t/handleVariantRecord.jl 2>&1 | tail -20
```

Expected: All tests pass.

- [ ] **Step 7: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "$(cat <<'EOF'
refactor: strip variationFeature.dat to 14 genomic columns

Move write_snp_feature outside annotation loop; drop transcript_id and
all CDS-specific columns (products, ref_codon, pos_in_cds, DFS IDs).
One row per genomic position instead of one per transcript.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Merge product.dat + cache.tsv into transcript_product.dat

**Files:**
- Modify: `testing/t/handleVariantRecord.jl` (update write_product_file tests)
- Modify: `bin/processSequenceVariations.jl` — `OutputWriters` struct, `open_output_writers`, `close_output_writers`, `write_product_file` → `write_transcript_product`, remove `write_cache_entries`, update call sites, update `main()`

### New column set for transcript_product.dat (12, 0-indexed)

| idx | name | note |
|---|---|---|
| 0 | seq_id | matches cache col 1 — `peek_sort_key` works unchanged |
| 1 | location | matches cache col 2 |
| 2 | transcript_id | matches cache col 3 |
| 3 | pos_in_cds | matches cache col 4 |
| 4 | pos_in_protein | `div(pos_in_cds - 1, 3) + 1` |
| 5 | codon | was col 2 |
| 6 | pos_in_codon | was col 3 |
| 7 | count | was col 5 |
| 8 | product | was col 6 |
| 9 | matches_ref_codon | was col 7 |
| 10 | matches_ref_product | was col 8 |
| 11 | downstream_of_frameshift_strain_ids | moved from variationFeature.dat |

Header starts with `#` so `open_cache_peeked` skips it (it skips `#`-prefixed lines).

- [ ] **Step 1: Update existing `write_product_file` unit tests to the new name and schema**

In `testing/t/handleVariantRecord.jl`, replace every occurrence of `write_product_file(buf, ...)` with `write_transcript_product(buf, ..., Dict{String,Int}())`. Update column-index assertions throughout the `write_product_file` test section:

```julia
# Old assertions used length(fields) == 9, fields[3], fields[6], fields[8], fields[9]
# New: length == 12; codon=fields[6], count=fields[8], matches_ref_codon=fields[10], matches_ref_product=fields[11]
```

Concretely, update the existing 6 `write_product_file` tests:

```julia
# "write_product_file skips DFS strains..."
buf = IOBuffer()
write_transcript_product(buf, [v1, v2], ann, "chr1", 100, Dict{String,Int}())
lines = filter(!isempty, split(String(take!(buf)), "\n"))
@test length(lines) == 1
fields = split(lines[1], "\t")
@test length(fields) == 12
@test fields[6] == "ATT"          # codon is now col 6 (1-indexed)

# "deduplicates identical codons..."
# fields[6] == "ATT"; fields[8] == "2"  (count)

# "matches_ref_codon=1..."
# fields[10] == "1"; fields[11] == "1"

# "matches_ref_codon=0, matches_ref_product=0..."
# fields[10] == "0"; fields[11] == "0"

# "synonymous..."
# fields[10] == "0"; fields[11] == "1"

# "has 9 columns per row" → change to 12
@test length(split(lines[1], "\t")) == 12
```

Also add tests for the new columns:
```julia
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
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
julia testing/t/handleVariantRecord.jl 2>&1 | tail -20
```

Expected: `UndefVarError: write_transcript_product not defined`

- [ ] **Step 3: Update `OutputWriters` struct (~line 422)**

```julia
struct OutputWriters
    vcf_fh::IO
    snp_fh::IO
    allele_fh::IO
    transcript_product_fh::IO
end
```

- [ ] **Step 4: Update `open_output_writers` (~line 1071–1085)**

```julia
function open_output_writers(output_vcf::String)
    vcf_fh = open(output_vcf, "w")
    snp_fh = open("snpFeature.dat", "w")
    write(snp_fh, "location\tseq_id\treference_strain\tref_allele\tmajor_allele\tminor_allele\tmajor_allele_strain_count\tminor_allele_strain_count\tmajor_allele_frequency\tminor_allele_frequency\tdistinct_strain_count\tdistinct_allele_count\ttotal_ploidy_count\tis_coding\n")
    allele_fh = open("allele.dat", "w")
    write(allele_fh, "location\tseq_id\tallele\tdistinct_strain_count\tallele_frequency\tavg_coverage\tavg_percent\tstrain_ids\tmatches_reference\n")
    tp_fh = open("transcript_product.dat", "w")
    write(tp_fh, "#seq_id\tlocation\ttranscript_id\tpos_in_cds\tpos_in_protein\tcodon\tpos_in_codon\tcount\tproduct\tmatches_ref_codon\tmatches_ref_product\tdownstream_of_frameshift_strain_ids\n")
    OutputWriters(vcf_fh, snp_fh, allele_fh, tp_fh)
end
```

- [ ] **Step 5: Update `close_output_writers` (~line 1092)**

```julia
function close_output_writers(writers::OutputWriters)
    close(writers.vcf_fh)
    close(writers.snp_fh)
    close(writers.allele_fh)
    close(writers.transcript_product_fh)
end
```

- [ ] **Step 6: Rename `write_product_file` → `write_transcript_product` and update signature/body**

Replace the function signature and the `write(...)` inside the codon loop:

```julia
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
```

- [ ] **Step 7: Delete `write_cache_entries` function (~line 719–724)**

Remove the entire function:
```julia
function write_cache_entries(fh::IO, chrom::String, pos::Int, annotations::Vector{PositionAnnotation})
    ...
end
```

- [ ] **Step 8: Update call sites in `handle_variant_record!`**

Remove:
```julia
write_cache_entries(writers.cache_fh, seq_id, location, annotations)
```

Replace:
```julia
write_product_file(writers.product_fh, all_vars, annotation, seq_id, location)
```
with:
```julia
write_transcript_product(writers.transcript_product_fh, all_vars, annotation, seq_id, location, ctx.sample_id_map)
```

- [ ] **Step 9: Update `main()` (~line 1943)**

```julia
writers = open_output_writers(args["output_vcf"])
```

(Remove the `args["output_cache"]` argument.)

- [ ] **Step 10: Run tests to verify they pass**

```bash
julia testing/t/handleVariantRecord.jl 2>&1 | tail -20
```

Expected: All tests pass.

- [ ] **Step 11: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "$(cat <<'EOF'
refactor: merge product.dat+cache.tsv into transcript_product.dat

New 12-column file: seq_id, location, transcript_id, pos_in_cds,
pos_in_protein, codon, pos_in_codon, count, product, matches_ref_codon,
matches_ref_product, downstream_of_frameshift_strain_ids. First 4 cols
preserve cache sort-key order so parse_cache_tsv_record works unchanged.
Removes write_cache_entries and the separate output_cache argument.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Add HSSS helpers and unit tests

**Files:**
- Modify: `testing/t/handleVariantRecord.jl` (add HSSS unit tests)
- Modify: `bin/processSequenceVariations.jl` (add constants, struct, open/close/write functions)

### Binary record: 8 bytes matching Perl `pack("slcc", seq_idx, location, allele_code, product_code)`
- `Int16` (little-endian) — seq_index
- `Int32` (little-endian) — location
- `Int8` — allele_code
- `Int8` — product_code

- [ ] **Step 1: Write failing unit tests for HSSS helpers**

Append to `testing/t/handleVariantRecord.jl`:

```julia
# ---------------------------------------------------------------------------
# HSSS binary strain files
# ---------------------------------------------------------------------------

@testset "open_hsss_writers creates directories and strainIdToName.dat" begin
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
    end
end

@testset "write_hsss_position writes binary record and contigIdToSourceId" begin
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1"], base)

    # ref: A, s1: T (passes all cutoffs)
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"
    v_ref.percent = "100"; v_ref.matches_reference = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "T"; v_s1.reference  = "A"
    v_s1.percent = "75"; v_s1.matches_reference = 0

    write_hsss_position!(state, [v_ref, v_s1], "ref", "LmjF.01", 200, ["s1"], Int8(77))
    close_hsss_writers(state)

    # Verify contig map written to all 4 dirs
    for cutoff in [20, 40, 60, 80]
        dir = joinpath(base, "hsss_readFreq$(cutoff)")
        contig_lines = readlines(joinpath(dir, "contigIdToSourceId.dat"))
        @test contig_lines[1] == "1\tLmjF.01"
    end

    # Verify binary record in readFreq20/referenceGenome.dat
    ref_bytes = read(joinpath(base, "hsss_readFreq20", "referenceGenome.dat"))
    @test length(ref_bytes) == 8
    @test ltoh(reinterpret(Int16, ref_bytes[1:2])[1]) == Int16(1)   # seq_index
    @test ltoh(reinterpret(Int32, ref_bytes[3:6])[1]) == Int32(200) # location
    @test reinterpret(Int8, ref_bytes[7:7])[1] == Int8(1)           # A = 1
    @test reinterpret(Int8, ref_bytes[8:8])[1] == Int8(77)          # product_code

    # Verify binary record for s1 (index 2)
    s1_bytes = read(joinpath(base, "hsss_readFreq20", "2"))
    @test length(s1_bytes) == 8
    @test reinterpret(Int8, s1_bytes[7:7])[1] == Int8(4)   # T = 4
    @test reinterpret(Int8, s1_bytes[8:8])[1] == Int8(77)
end

@testset "write_hsss_position skips position if only ref-matching strains at cutoff" begin
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1"], base)

    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"
    v_ref.percent = "100"; v_ref.matches_reference = 1
    # s1 matches reference
    v_s1 = Variation(); v_s1.strain = "s1"; v_s1.base = "A"; v_s1.reference = "A"
    v_s1.percent = "100"; v_s1.matches_reference = 1

    write_hsss_position!(state, [v_ref, v_s1], "ref", "chr1", 100, ["s1"], Int8(77))
    close_hsss_writers(state)

    # No records written anywhere
    @test filesize(joinpath(base, "hsss_readFreq20", "referenceGenome.dat")) == 0
    @test filesize(joinpath(base, "hsss_readFreq20", "2")) == 0
end

@testset "write_hsss_position writes unknown record for absent strain" begin
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1", "s2"], base)

    # Only ref and s1 present; s2 is absent
    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"
    v_ref.percent = "100"; v_ref.matches_reference = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "T"; v_s1.reference  = "A"
    v_s1.percent = "60"; v_s1.matches_reference = 0

    write_hsss_position!(state, [v_ref, v_s1], "ref", "chr1", 50, ["s1", "s2"], Int8(0))
    close_hsss_writers(state)

    # s2 (index 3) gets unknown record
    s2_bytes = read(joinpath(base, "hsss_readFreq20", "3"))
    @test length(s2_bytes) == 8
    @test reinterpret(Int8, s2_bytes[7:7])[1] == Int8(0)  # unknown allele
    @test reinterpret(Int8, s2_bytes[8:8])[1] == Int8(0)  # unknown product
end

@testset "write_hsss_position filters by cutoff: s1 at 35% skips readFreq40" begin
    base = mktempdir()
    state = open_hsss_writers("ref", ["s1"], base)

    v_ref = Variation(); v_ref.strain = "ref"; v_ref.base = "A"; v_ref.reference = "A"
    v_ref.percent = "100"; v_ref.matches_reference = 1
    v_s1  = Variation(); v_s1.strain  = "s1";  v_s1.base  = "T"; v_s1.reference  = "A"
    v_s1.percent = "35"; v_s1.matches_reference = 0   # passes 20, fails 40

    write_hsss_position!(state, [v_ref, v_s1], "ref", "chr1", 100, ["s1"], Int8(0))
    close_hsss_writers(state)

    # readFreq20: s1 passes → record written
    @test filesize(joinpath(base, "hsss_readFreq20", "2")) == 8
    # readFreq40: s1 fails → position skipped (no record, not even unknown, since no non-ref passes)
    @test filesize(joinpath(base, "hsss_readFreq40", "2")) == 0
end
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
julia testing/t/handleVariantRecord.jl 2>&1 | tail -20
```

Expected: `UndefVarError: open_hsss_writers not defined`

- [ ] **Step 3: Add constants and `HsssState` struct to `processSequenceVariations.jl`**

Add after the `IUPAC_COMPRESS` constant block (~line 88):

```julia
# ---------------------------------------------------------------------------
# HSSS (HighSpeedSnpSearch) constants and state
# ---------------------------------------------------------------------------

const HSSS_CUTOFFS = [20, 40, 60, 80]

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
```

- [ ] **Step 4: Add `open_hsss_writers` function**

Add after `HsssState` definition:

```julia
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
```

- [ ] **Step 5: Add `close_hsss_writers` function**

```julia
function close_hsss_writers(state::HsssState)
    for fh in state.ref_fhs;    close(fh); end
    for fh in state.contig_fhs; close(fh); end
    for sfhs in state.strain_fhs
        for (_, fh) in sfhs; close(fh); end
    end
end
```

- [ ] **Step 6: Add `write_hsss_position!` function**

```julia
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

            passes = any(v -> !isempty(v.percent) && parse(Float64, v.percent) >= cutoff, svars)
            if !passes
                # Present but below cutoff: treat as unknown
                write_hsss_record(sfh, seq_idx, location, Int8(0), Int8(0))
                continue
            end

            is_het = length(svars) > 1
            for sv in svars
                # Skip non-het, reference-matching variations
                !is_het && sv.matches_reference == 1 && continue
                allele_c = get(HSSS_ALLELE_CODE, isempty(sv.base) ? ' ' : sv.base[1], Int8(0))
                write_hsss_record(sfh, seq_idx, location, allele_c, product_code)
            end
        end
    end
end

function write_hsss_record(fh::IO, seq_idx::Int, location::Int, allele_c::Int8, product_c::Int8)
    write(fh, htol(Int16(seq_idx)))
    write(fh, htol(Int32(location)))
    write(fh, allele_c)
    write(fh, product_c)
end
```

- [ ] **Step 7: Run tests to verify they pass**

```bash
julia testing/t/handleVariantRecord.jl 2>&1 | tail -20
```

Expected: All tests pass.

- [ ] **Step 8: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "$(cat <<'EOF'
feat: add HSSS binary strain file helpers

Add HsssState, open/close_hsss_writers, write_hsss_position!, and
write_hsss_record. Produces 4 readFreqXX directories per run, each
with binary strain files matching the HighSpeedSnpSearch format
(8 bytes/record: Int16 seq_idx, Int32 location, Int8 allele, Int8 product).

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Integrate HSSS into main loop and wire Nextflow

**Files:**
- Modify: `bin/processSequenceVariations.jl` — `OutputWriters`, `open_output_writers`, `close_output_writers`, `handle_variant_record!`, `main()`
- Modify: `modules/mergeExperiments.nf`

- [ ] **Step 1: Add `hsss` field to `OutputWriters` struct (~line 422)**

```julia
struct OutputWriters
    vcf_fh::IO
    snp_fh::IO
    allele_fh::IO
    transcript_product_fh::IO
    hsss::HsssState
end
```

- [ ] **Step 2: Update `open_output_writers` to accept and initialize HSSS**

```julia
function open_output_writers(output_vcf::String, reference_strain::String,
                              all_strains::Vector{String})
    vcf_fh = open(output_vcf, "w")
    snp_fh = open("snpFeature.dat", "w")
    write(snp_fh, "location\tseq_id\treference_strain\tref_allele\tmajor_allele\tminor_allele\tmajor_allele_strain_count\tminor_allele_strain_count\tmajor_allele_frequency\tminor_allele_frequency\tdistinct_strain_count\tdistinct_allele_count\ttotal_ploidy_count\tis_coding\n")
    allele_fh = open("allele.dat", "w")
    write(allele_fh, "location\tseq_id\tallele\tdistinct_strain_count\tallele_frequency\tavg_coverage\tavg_percent\tstrain_ids\tmatches_reference\n")
    tp_fh = open("transcript_product.dat", "w")
    write(tp_fh, "#seq_id\tlocation\ttranscript_id\tpos_in_cds\tpos_in_protein\tcodon\tpos_in_codon\tcount\tproduct\tmatches_ref_codon\tmatches_ref_product\tdownstream_of_frameshift_strain_ids\n")
    hsss = open_hsss_writers(reference_strain, all_strains)
    OutputWriters(vcf_fh, snp_fh, allele_fh, tp_fh, hsss)
end
```

- [ ] **Step 3: Update `close_output_writers`**

```julia
function close_output_writers(writers::OutputWriters)
    close(writers.vcf_fh)
    close(writers.snp_fh)
    close(writers.allele_fh)
    close(writers.transcript_product_fh)
    close_hsss_writers(writers.hsss)
end
```

- [ ] **Step 4: Integrate HSSS into `handle_variant_record!`**

Add before the annotation loop (alongside `first_all_vars = nothing`):
```julia
all_annotation_products = String[]
```

Inside the annotation loop, after `any_output = true`, add:
```julia
for v in variations
    append!(all_annotation_products, v.product)
end
```

After `write_snp_feature(...)` (the call added in Task 1), add:
```julia
unique_prods   = unique(all_annotation_products)
hsss_prod_code = length(unique_prods) == 1 ?
    Int8(codepoint(only(unique_prods)[1])) : Int8(0)
write_hsss_position!(writers.hsss, first_all_vars, ctx.reference_strain,
                     seq_id, location, all_strains, hsss_prod_code)
```

- [ ] **Step 5: Update `main()` call to `open_output_writers`**

```julia
writers = open_output_writers(args["output_vcf"], args["reference_strain"], all_strains)
```

- [ ] **Step 6: Update `modules/mergeExperiments.nf`**

In the `processSeqVars` process:

**Remove** these `publishDir` lines:
```nextflow
publishDir "$params.outputDir", mode: "copy", pattern: 'output_cache.tsv', saveAs: { 'cache.tsv' }
publishDir "$params.outputDir", mode: "copy", pattern: 'product.dat'
```

**Add/replace** with:
```nextflow
publishDir "$params.outputDir", mode: "copy", pattern: 'transcript_product.dat'
publishDir "$params.outputDir", mode: "copy", pattern: 'hsss_readFreq20'
publishDir "$params.outputDir", mode: "copy", pattern: 'hsss_readFreq40'
publishDir "$params.outputDir", mode: "copy", pattern: 'hsss_readFreq60'
publishDir "$params.outputDir", mode: "copy", pattern: 'hsss_readFreq80'
```

**In the `output:` block**, remove:
```nextflow
path 'output_cache.tsv', emit: outputCache
path 'product.dat',      emit: productFile
```

Add:
```nextflow
path 'transcript_product.dat', emit: transcriptProductFile
path 'hsss_readFreq20',        emit: hsssReadFreq20
path 'hsss_readFreq40',        emit: hsssReadFreq40
path 'hsss_readFreq60',        emit: hsssReadFreq60
path 'hsss_readFreq80',        emit: hsssReadFreq80
```

**In the `script:` block**, remove from the Julia CLI call:
```
      --output_cache output_cache.tsv \\
```

**In the `stub:` block**, remove:
```nextflow
touch output_cache.tsv
touch product.dat
```

Add:
```nextflow
touch transcript_product.dat
mkdir hsss_readFreq20
mkdir hsss_readFreq40
mkdir hsss_readFreq60
mkdir hsss_readFreq80
```

- [ ] **Step 7: Run Julia unit tests to verify nothing broke**

```bash
julia testing/t/handleVariantRecord.jl 2>&1 | tail -20
```

Expected: All tests pass.

- [ ] **Step 8: Commit**

```bash
git add bin/processSequenceVariations.jl modules/mergeExperiments.nf
git commit -m "$(cat <<'EOF'
feat: integrate HSSS into main loop and wire Nextflow outputs

Wire HsssState into OutputWriters; call write_hsss_position! after each
variant position with product_code derived from all annotations.
Remove --output_cache arg and cache.tsv output; add transcript_product.dat
and 4 hsss_readFreqXX directory outputs to processSeqVars process.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Update e2e tests

**Files:**
- Modify: `testing/t/test_mergeExperiments_e2e.py`

### variationFeature.dat — column index mapping (0-indexed)

Old → New: location[0]→[0], seq_id[2]→[1], reference_strain[3]→[2], ref_allele[4]→[3], major_allele[6]→[4], minor_allele[7]→[5], major_allele_strain_count[8]→[6], minor_allele_strain_count[9]→[7], major_allele_frequency[20]→[8], minor_allele_frequency[21]→[9], distinct_strain_count[12]→[10], distinct_allele_count[13]→[11], total_ploidy_count[15]→[12], is_coding[14]→[13].

### product.dat → transcript_product.dat — column index mapping (0-indexed)

Old → New: location[0]→[1] (now location is col 2, seq_id is col 1), seq_id[1]→[0], transcript_id[4]→[2], codon[2]→[5], pos_in_codon[3]→[6], count[5]→[7], product[6]→[8], matches_ref_codon[7]→[9], matches_ref_product[8]→[10].

- [ ] **Step 1: Update `_read_variation_feature` helper — no change needed** (still reads variationFeature.dat, still skips header line)

- [ ] **Step 2: Update variationFeature.dat tests**

```python
# test_variation_feature_column_count: 22 → 14
bad = [i + 1 for i, r in enumerate(rows) if len(r) != 14]
assert not bad, f"Rows with wrong column count (expected 14): {bad[:5]}"

# test_variation_feature_seq_id_nonempty: r[2] → r[1]
bad = [i + 1 for i, r in enumerate(rows) if not r[1].strip()]
assert not bad, f"Rows with empty seq_id (col 2): {bad[:5]}"

# test_variation_feature_reference_strain: r[3] → r[2]
values = {r[2] for r in rows}

# test_variation_feature_major_allele_nonempty: r[6] → r[4]
bad = [i + 1 for i, r in enumerate(rows) if not r[4].strip()]

# test_variation_feature_major_allele_strain_count_positive: r[8] → r[6]
bad = [i + 1 for i, r in enumerate(rows) if not r[6].isdigit() or int(r[6]) <= 0]

# test_variation_feature_distinct_strain_count_in_range: r[12] → r[10]
bad = [i + 1 for i, r in enumerate(rows)
       if not r[10].isdigit() or not (1 <= int(r[10]) <= n_strains + 1)]

# test_variation_feature_is_coding_binary: r[14] → r[13]
bad = [i + 1 for i, r in enumerate(rows) if r[13] not in ('0', '1')]

# test_variation_feature_ref_allele_nonempty: r[4] → r[3]
bad = [i + 1 for i, r in enumerate(rows) if not r[3].strip()]

# test_variation_feature_distinct_allele_count_gte_1: r[13] → r[11]
bad = [i + 1 for i, r in enumerate(rows) if not r[11].isdigit() or int(r[11]) < 1]

# test_variation_feature_total_ploidy_count_gte_strain_count: r[15] → r[12], r[12] → r[10]
bad = [i + 1 for i, r in enumerate(rows)
       if not r[12].isdigit() or int(r[12]) < int(r[10])]

# test_variation_feature_major_allele_frequency_in_range: r[20] → r[8]
f = float(r[8])

# test_variation_feature_allele_frequencies_sum_to_one: r[20] → r[8], r[21] → r[9]
major = float(r[8])
minor = float(r[9]) if r[9] != '' else 0.0

# test_variation_feature_no_iupac_alleles: r[6], r[7] → r[4], r[5]
bad = [(i + 1, r[4], r[5]) for i, r in enumerate(rows)
       if iupac_pattern.match(r[4]) or (r[5] and iupac_pattern.match(r[5]))]

# test_reference_strain_nonempty_and_consistent_with_vf: r[3] → r[2]
ref_strain = vf_rows[0][2]
inconsistent = [i + 1 for i, r in enumerate(vf_rows) if r[2] != ref_strain]
```

- [ ] **Step 3: Delete tests that reference dropped columns**

Remove these test functions entirely:
- `test_variation_feature_has_nonsynonymous_binary` (col 5 dropped)
- `test_variation_feature_coding_rows_have_transcript` (transcript_id dropped)
- `test_variation_feature_has_stop_codon_binary` (col 16 dropped)
- `test_variation_feature_nonsynonymous_implies_coding` (has_nonsynonymous dropped)
- `test_variation_feature_dfs_strain_ids_format` (DFS col dropped from this file)
- `test_variation_feature_pos_in_cds_coding_only` (pos_in_cds dropped)

- [ ] **Step 4: Update product.dat → transcript_product.dat**

Rename the section, reader helper, and all test functions from `product` to `transcript_product`:

```python
def _read_transcript_product(work_dirs):
    path = os.path.join(work_dirs['processSeqVars'], 'transcript_product.dat')
    rows = []
    with open(path) as f:
        next(f)  # skip header (starts with #)
        for line in f:
            rows.append(line.rstrip('\n').split('\t'))
    return rows

def test_transcript_product_dat_exists(work_dirs):
    assert os.path.exists(
        os.path.join(work_dirs['processSeqVars'], 'transcript_product.dat'))

def test_transcript_product_dat_column_count(work_dirs):
    rows = _read_transcript_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if len(r) != 12]
    assert not bad, f"Rows with wrong column count (expected 12): {bad[:5]}"

def test_transcript_product_dat_location_positive_int(work_dirs):
    rows = _read_transcript_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[1].isdigit() or int(r[1]) <= 0]
    assert not bad, f"Rows with invalid location (col 2): {bad[:5]}"

def test_transcript_product_dat_seq_id_nonempty(work_dirs):
    rows = _read_transcript_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[0].strip()]
    assert not bad, f"Rows with empty seq_id (col 1): {bad[:5]}"

def test_transcript_product_dat_transcript_id_nonempty(work_dirs):
    rows = _read_transcript_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[2].strip()]
    assert not bad, f"Rows with empty transcript_id (col 3): {bad[:5]}"

def test_transcript_product_dat_pos_in_cds_positive(work_dirs):
    rows = _read_transcript_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[3].isdigit() or int(r[3]) <= 0]
    assert not bad, f"Rows with pos_in_cds <= 0 (col 4): {bad[:5]}"

def test_transcript_product_dat_pos_in_protein_positive(work_dirs):
    rows = _read_transcript_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[4].isdigit() or int(r[4]) <= 0]
    assert not bad, f"Rows with pos_in_protein <= 0 (col 5): {bad[:5]}"

def test_transcript_product_dat_codon_is_three_acgt(work_dirs):
    pattern = re.compile(r'^[ACGTN]{3}$')
    rows = _read_transcript_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not pattern.match(r[5])]
    assert not bad, f"Rows with invalid codon (col 6): {bad[:5]}"

def test_transcript_product_dat_pos_in_codon_1_2_3(work_dirs):
    rows = _read_transcript_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if r[6] not in ('1', '2', '3')]
    assert not bad, f"Rows with pos_in_codon not in {{1,2,3}} (col 7): {bad[:5]}"

def test_transcript_product_dat_count_non_negative(work_dirs):
    rows = _read_transcript_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[7].lstrip('-').isdigit() or int(r[7]) < 0]
    assert not bad, f"Rows with count < 0 (col 8): {bad[:5]}"

def test_transcript_product_dat_amino_acid_single_char_or_stop(work_dirs):
    rows = _read_transcript_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if len(r[8]) != 1]
    assert not bad, f"Rows where amino_acid not single char (col 9): {bad[:5]}"

def test_transcript_product_dat_no_duplicate_rows(work_dirs):
    rows = _read_transcript_product(work_dirs)
    seen = set()
    dups = []
    for i, r in enumerate(rows):
        key = tuple(r)
        if key in seen:
            dups.append(i + 1)
        seen.add(key)
    assert not dups, f"Duplicate rows in transcript_product.dat: {dups[:5]}"
```

- [ ] **Step 5: Update `test_product_dat_transcripts_in_coding_sequences` to use new reader**

```python
def test_transcript_product_transcripts_in_coding_sequences(work_dirs):
    db_path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    with sqlite3.connect(db_path) as conn:
        known = {row[0] for row in conn.execute(
            "SELECT DISTINCT transcript_id FROM coding_sequences")}
    rows = _read_transcript_product(work_dirs)
    transcripts = {r[2] for r in rows if r[2]}   # col 2 = transcript_id
    orphans = transcripts - known
    assert not orphans, \
        f"transcript_product.dat references unknown transcripts: {list(orphans)[:5]}"
```

- [ ] **Step 6: Add HSSS e2e tests**

```python
# ---------------------------------------------------------------------------
# Layer 1: processSeqVars — HSSS binary strain files
# ---------------------------------------------------------------------------

def test_hsss_directories_exist(work_dirs):
    for cutoff in [20, 40, 60, 80]:
        path = os.path.join(work_dirs['processSeqVars'], f'hsss_readFreq{cutoff}')
        assert os.path.isdir(path), f"hsss_readFreq{cutoff} directory missing"


def test_hsss_strain_map_files_exist(work_dirs):
    for cutoff in [20, 40, 60, 80]:
        d = os.path.join(work_dirs['processSeqVars'], f'hsss_readFreq{cutoff}')
        assert os.path.exists(os.path.join(d, 'strainIdToName.dat'))
        assert os.path.exists(os.path.join(d, 'contigIdToSourceId.dat'))
        assert os.path.exists(os.path.join(d, 'referenceGenome.dat'))


def test_hsss_strain_map_consistent_across_cutoffs(work_dirs):
    """strainIdToName.dat must be identical across all 4 cutoffs."""
    contents = []
    for cutoff in [20, 40, 60, 80]:
        p = os.path.join(work_dirs['processSeqVars'],
                         f'hsss_readFreq{cutoff}', 'strainIdToName.dat')
        with open(p) as f:
            contents.append(f.read())
    assert len(set(contents)) == 1, "strainIdToName.dat differs across cutoffs"


def test_hsss_ref_genome_size_multiple_of_8(work_dirs):
    """referenceGenome.dat must be a multiple of 8 bytes (one record each)."""
    for cutoff in [20, 40, 60, 80]:
        p = os.path.join(work_dirs['processSeqVars'],
                         f'hsss_readFreq{cutoff}', 'referenceGenome.dat')
        size = os.path.getsize(p)
        assert size % 8 == 0, \
            f"referenceGenome.dat size {size} not multiple of 8 (readFreq{cutoff})"


def test_hsss_strain_files_size_multiple_of_8(work_dirs):
    """All per-strain binary files must be multiples of 8 bytes."""
    import struct
    for cutoff in [20, 40, 60, 80]:
        d = os.path.join(work_dirs['processSeqVars'], f'hsss_readFreq{cutoff}')
        strain_map = os.path.join(d, 'strainIdToName.dat')
        with open(strain_map) as f:
            indices = [line.strip().split('\t')[0] for line in f if line.strip()]
        for idx in indices[1:]:   # skip ref (index 1)
            p = os.path.join(d, idx)
            if os.path.exists(p):
                size = os.path.getsize(p)
                assert size % 8 == 0, \
                    f"Strain file {idx} size {size} not multiple of 8 (readFreq{cutoff})"


def test_hsss_higher_cutoff_no_more_records(work_dirs):
    """readFreq40 referenceGenome.dat must have <= records than readFreq20."""
    def record_count(path):
        return os.path.getsize(path) // 8 if os.path.exists(path) else 0

    d = work_dirs['processSeqVars']
    r20 = record_count(os.path.join(d, 'hsss_readFreq20', 'referenceGenome.dat'))
    r40 = record_count(os.path.join(d, 'hsss_readFreq40', 'referenceGenome.dat'))
    assert r40 <= r20, f"readFreq40 ({r40}) has more records than readFreq20 ({r20})"
```

- [ ] **Step 7: Commit**

```bash
git add testing/t/test_mergeExperiments_e2e.py
git commit -m "$(cat <<'EOF'
test: update e2e tests for variationFeature/transcript_product/HSSS

variationFeature.dat: 22 → 14 columns; update all column indices, remove
tests for dropped columns. product.dat → transcript_product.dat: rename
reader/tests, update column indices for 12-col format. Add HSSS directory
existence, file size, and cross-cutoff consistency checks.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Spec Coverage Self-Check

| Spec requirement | Task |
|---|---|
| variationFeature.dat — 14 genomic columns, one row per position | Task 1 |
| variationFeature.dat — drop transcript_id, products, ref_codon, pos_in_cds, DFS | Task 1 |
| transcript_product.dat — replaces product.dat + cache.tsv | Task 2 |
| transcript_product.dat — seq_id first (cache sort-key compatible) | Task 2 |
| transcript_product.dat — pos_in_protein added | Task 2 |
| transcript_product.dat — downstream_of_frameshift_strain_ids | Task 2 |
| transcript_product.dat — header starts with `#` | Task 2 |
| Remove `--output_cache` CLI arg | Task 2 |
| HSSS_ALLELE_CODE constant, HsssState struct | Task 3 |
| 4 cutoffs, 4 directories, correct file structure | Task 3 |
| Binary format: Int16/Int32/Int8/Int8, little-endian | Task 3 |
| Product encoding: unique across all annotations or 0 | Task 4 |
| write_hsss_position! in main loop | Task 4 |
| Nextflow: remove cache output, rename product, add 4 HSSS dirs | Task 4 |
| E2e test updates for all three file changes | Task 5 |
