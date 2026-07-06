# variationFeature per-class schema + allele identity fix — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give `variationFeature.dat` a per-class (SNP/indel) 31-column schema and re-key all locus allele aggregation by `(ref_span, allele)` instead of bare base string, fixing the deferred `matches_reference` bug and Finding 2 (deletion collapses into reference) at their shared root.

**Architecture:** Introduce one pure aggregation helper keyed by `(reference, base)` tuples that both `write_snp_feature` and `write_allele_file` consume (removing today's duplicated, divergent counting). `write_snp_feature` partitions keys into reference/SNP/indel buckets and emits per-class columns; `write_allele_file` emits one row per `(ref,base)` with per-row-correct `matches_reference`. Genomic HGVS is computed directly from each key's own ref span, retiring the collision-prone `build_allele_refs`/`allele_genomic_hgvs`.

**Tech Stack:** Julia 1.10 (`Test` stdlib), Python 3 (`pytest`) for e2e.

**Spec:** `docs/superpowers/specs/2026-07-06-variationfeature-per-class-schema-design.md`

## How to run the Julia tests (julia is NOT on PATH)
```
SCR="/tmp/claude-1000/-home-jbrestel-workspaces-dataLoad-project-home-dnaseq-nextflow/9391ce2c-fcc8-4210-8c04-1b0a8911e723/scratchpad"
JULIA_DEPOT_PATH="$SCR/julia_depot" "$SCR/julia-1.10.10/bin/julia" testing/t/handleVariantRecord.jl
```

## File structure

| File | Responsibility | Action |
|---|---|---|
| `bin/processSequenceVariations.jl` | aggregation helper, `write_snp_feature`, `write_allele_file`, header | Modify |
| `testing/t/handleVariantRecord.jl` | Julia unit tests | Modify (add tests) |
| `testing/t/test_mergeExperiments_e2e.py` | e2e column/invariant assertions | Modify |

## Key facts (verified against current source)
- `Variation` carries `reference`, `base`, `alt_allele` (het), `ploidy`, `coverage`, `percent`, `strain`. Het calls have non-empty `alt_allele`; hom calls have empty `alt_allele` and use `base`.
- Today `compute_allele_weight_map` (line ~1447) keys weights by bare `base`/`alt_allele`; `write_snp_feature` (~1463) and `write_allele_file` (~1567) key by bare string and take `ref_allele = variations[1].reference`. This is the collapse.
- `genomic_hgvs(seq_id, pos, ref, alt)` (~1823) already produces correct g. HGVS for any (ref,alt) — call it directly per key.
- Header for `variationFeature.dat` is written at line ~1255; `allele.dat` header at ~1257.

---

## Task 1: `(ref,base)`-keyed aggregation helper + classifier

**Files:**
- Modify: `bin/processSequenceVariations.jl` (add `AlleleStat`, `aggregate_locus_alleles`, `classify_allele` near `compute_allele_weight_map`, ~line 1447)
- Test: `testing/t/handleVariantRecord.jl`

- [ ] **Step 1: Write the failing test.** Append to `testing/t/handleVariantRecord.jl`:

```julia
# ---------------------------------------------------------------------------
# aggregate_locus_alleles / classify_allele — (ref,base)-keyed aggregation
# ---------------------------------------------------------------------------

function mkvar(; strain, reference, base, alt_allele="", ploidy=1,
                coverage="30", percent="100", matches_reference=0)
    v = Variation()
    v.strain = strain; v.reference = reference; v.base = base
    v.alt_allele = alt_allele; v.ploidy = ploidy
    v.coverage = coverage; v.percent = percent
    v.matches_reference = matches_reference
    v
end

@testset "classify_allele distinguishes reference/snp/indel" begin
    @test classify_allele("A", "A")     == :reference
    @test classify_allele("A", "G")     == :snp
    @test classify_allele("ACA", "A")   == :indel   # deletion
    @test classify_allele("A", "AT")    == :indel   # insertion
end

@testset "aggregate_locus_alleles keys by (ref,base), no collapse" begin
    # Mixed locus: SNP A>G (S1), deletion ACA>A (S2), reference genome (REF)
    vars = [
        mkvar(strain="S1", reference="A",   base="G"),
        mkvar(strain="S2", reference="ACA", base="A"),
        mkvar(strain="REF", reference="A",  base="A", matches_reference=1),
    ]
    (stats, total) = aggregate_locus_alleles(vars)
    @test total == 3
    # deletion A (ref ACA) and reference A (ref A) are DISTINCT keys
    @test haskey(stats, ("ACA", "A"))
    @test haskey(stats, ("A", "A"))
    @test haskey(stats, ("A", "G"))
    @test stats[("ACA","A")].weight == 1
    @test length(stats[("A","G")].strains) == 1
end
```

- [ ] **Step 2: Run to verify it fails.** Run the julia command above. Expected: FAIL — `classify_allele`/`aggregate_locus_alleles` not defined.

- [ ] **Step 3: Implement.** In `bin/processSequenceVariations.jl`, immediately before `function compute_allele_weight_map` (~line 1447), add:

```julia
mutable struct AlleleStat
    weight::Int
    strains::Set{String}
    cov_sum::Float64
    pct_sum::Float64
    entry_count::Int
end
AlleleStat() = AlleleStat(0, Set{String}(), 0.0, 0.0, 0)

"""
    classify_allele(ref, allele) -> Symbol

:reference if allele == ref; :snp if same length (substitution); :indel otherwise.
"""
function classify_allele(ref::String, allele::String)::Symbol
    allele == ref                 ? :reference :
    length(allele) == length(ref) ? :snp : :indel
end

"""
    aggregate_locus_alleles(variations) -> (Dict{Tuple{String,String},AlleleStat}, total_weight)

Ploidy-weighted allele aggregation keyed by the (reference_span, allele) tuple, so
a deletion product (e.g. "A" from ref "ACA") never collides with a same-string SNP
reference. Het calls split ploidy across their reference and alt components exactly
as write_allele_file did. Shared by write_snp_feature and write_allele_file.
"""
function aggregate_locus_alleles(variations::Vector{Variation})::Tuple{Dict{Tuple{String,String},AlleleStat}, Int}
    stats = Dict{Tuple{String,String},AlleleStat}()
    total = 0
    add! = function(ref, allele, weight, strain, cov, pct)
        st = get!(stats, (ref, allele), AlleleStat())
        st.weight      += weight
        push!(st.strains, strain)
        st.cov_sum     += cov
        st.pct_sum     += pct
        st.entry_count += 1
    end
    for v in variations
        cov = isempty(v.coverage) ? 0.0 : parse(Float64, v.coverage)
        pct = isempty(v.percent)  ? 0.0 : parse(Float64, v.percent)
        if !isempty(v.alt_allele)                    # het: ref + alt components
            add!(v.reference, v.reference, 1, v.strain, cov, 100.0 - pct)
            add!(v.reference, v.alt_allele, 1, v.strain, cov, pct)
            total += 2
        else                                         # hom
            add!(v.reference, v.base, v.ploidy, v.strain, cov, pct)
            total += v.ploidy
        end
    end
    (stats, total)
end
```

- [ ] **Step 4: Run to verify it passes.** Run the julia command. Expected: PASS (new testsets + no regressions).

- [ ] **Step 5: Commit.**
```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "$(printf 'annotator: add (ref,base)-keyed locus allele aggregation helper\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 2: Rewrite `write_snp_feature` to the 31-column per-class schema

**Files:**
- Modify: `bin/processSequenceVariations.jl` — `write_snp_feature` (~1463–1565) and the header write (~line 1255)
- Test: `testing/t/handleVariantRecord.jl`

- [ ] **Step 1: Write the failing test.** Append to `testing/t/handleVariantRecord.jl`:

```julia
@testset "write_snp_feature emits per-class SNP+indel columns without collapse" begin
    # Mixed locus: SNP A>G (S1), deletion ACA>A (S2), reference genome (REF)
    vars = [
        mkvar(strain="S1", reference="A",   base="G", coverage="30", percent="100"),
        mkvar(strain="S2", reference="ACA", base="A", coverage="20", percent="100"),
        mkvar(strain="REF", reference="A",  base="A", coverage="0",  percent="100", matches_reference=1),
    ]
    buf = IOBuffer()
    write_snp_feature(buf, vars, 1, "LmjF.01", 13850, "REF", ["S1","S2"])
    cols = split(strip(String(take!(buf))), '\t')
    @test length(cols) == 31
    # locus-invariant
    @test cols[1] == "13850"          # location
    @test cols[5] == "MIXED"          # variant_type
    # SNP family (cols 13-21): ref A, major G
    @test cols[13] == "A"             # snp_ref_allele
    @test cols[14] == "G"             # snp_major_allele
    @test cols[20] == "LmjF.01:g.13850A>G"   # snp_major_genomic_hgvs
    # indel family (cols 22-31): ref ACA, major A (deletion), frame effect
    @test cols[22] == "ACA"           # indel_ref_allele
    @test cols[23] == "A"             # indel_major_allele
    @test occursin("del", cols[30])   # indel_major_genomic_hgvs
    @test cols[31] == "frameshift"    # indel_frame_effect (len 1-3=-2)
    # the deletion did NOT collapse into the SNP family or reference
    @test cols[14] != "A"
end

@testset "write_snp_feature SNP-only locus leaves indel family empty" begin
    vars = [
        mkvar(strain="S1", reference="A", base="G"),
        mkvar(strain="REF", reference="A", base="A", matches_reference=1),
    ]
    buf = IOBuffer()
    write_snp_feature(buf, vars, 1, "chr1", 100, "REF", ["S1"])
    cols = split(strip(String(take!(buf))), '\t')
    @test cols[5] == "SNV"
    @test cols[22] == ""              # indel_ref_allele empty
    @test cols[31] == ""              # indel_frame_effect empty
end
```

- [ ] **Step 2: Run to verify it fails.** Run the julia command. Expected: FAIL — old `write_snp_feature` emits 23 columns in the old order.

- [ ] **Step 3: Replace the header.** In `bin/processSequenceVariations.jl` at ~line 1255, replace the `variationFeature` header `write(snp_fh, "location\t...\n")` line with:

```julia
    write(snp_fh, "location\tseq_id\treference_strain\tis_coding\tvariant_type\tdistinct_strain_count\tcalled_strain_count\tno_call_strain_count\tcall_rate\ttotal_ploidy_count\tref_allele_frequency\thet_strain_count\tsnp_ref_allele\tsnp_major_allele\tsnp_major_allele_frequency\tsnp_major_allele_strain_count\tsnp_minor_allele\tsnp_minor_allele_frequency\tsnp_minor_allele_strain_count\tsnp_major_genomic_hgvs\tsnp_minor_genomic_hgvs\tindel_ref_allele\tindel_major_allele\tindel_major_allele_frequency\tindel_major_allele_strain_count\tindel_minor_allele\tindel_minor_allele_frequency\tindel_minor_allele_strain_count\tindel_major_genomic_hgvs\tindel_minor_genomic_hgvs\tindel_frame_effect\n")
```

- [ ] **Step 4: Replace `write_snp_feature`.** Replace the entire function body (~1463–1565) with:

```julia
function write_snp_feature(
    snp_fh::IO,
    variations::Vector{Variation},
    is_coding::Int,
    seq_id::String,
    location::Int,
    reference_strain::String,
    sequenced_strains::Vector{String}
)
    (stats, total_weight) = aggregate_locus_alleles(variations)

    ref_weight = 0
    snp_keys   = Tuple{String,String}[]
    indel_keys = Tuple{String,String}[]
    for (key, st) in stats
        c = classify_allele(key[1], key[2])
        if c == :reference
            ref_weight += st.weight
        elseif c == :snp
            push!(snp_keys, key)
        else
            push!(indel_keys, key)
        end
    end
    rank(ks) = sort(ks; by = k -> (-stats[k].weight, k[2]))
    snp_keys   = rank(snp_keys)
    indel_keys = rank(indel_keys)

    # locus-invariant
    strain_set = Set{String}(v.strain for v in variations)
    distinct_strain_count = length(strain_set)
    called_strain_count   = count(s -> s != reference_strain, strain_set)
    total_sequenced       = length(sequenced_strains)
    no_call_strain_count  = max(0, total_sequenced - called_strain_count)
    call_rate = total_sequenced > 0 ?
        @sprintf("%.4f", called_strain_count / total_sequenced) : ""
    het_strain_count = count(v -> !isempty(v.alt_allele), variations)
    has_snp   = !isempty(snp_keys)
    has_indel = !isempty(indel_keys)
    variant_type = has_snp && has_indel ? "MIXED" : (has_indel ? "INDEL" : "SNV")
    ref_allele_frequency = @sprintf("%.4f", ref_weight / total_weight)

    # per-class fields: (ref, major, major_freq, major_cnt, minor, minor_freq,
    # minor_cnt, major_hgvs, minor_hgvs) — all "" when the class is absent
    class_fields = function(keys)
        isempty(keys) && return ("","","","","","","","","")
        ref  = keys[1][1]
        maj  = keys[1][2]
        majf = @sprintf("%.4f", stats[keys[1]].weight / total_weight)
        majc = string(length(stats[keys[1]].strains))
        majh = genomic_hgvs(seq_id, location, ref, maj)
        if length(keys) > 1
            mn   = keys[2][2]
            mnf  = @sprintf("%.4f", stats[keys[2]].weight / total_weight)
            mnc  = string(length(stats[keys[2]].strains))
            mnh  = genomic_hgvs(seq_id, location, keys[2][1], mn)
        else
            mn = ""; mnf = ""; mnc = ""; mnh = ""
        end
        (ref, maj, majf, majc, mn, mnf, mnc, majh, mnh)
    end
    (sr, smaj, smajf, smajc, smin, sminf, sminc, smajh, sminh) = class_fields(snp_keys)
    (ir, imaj, imajf, imajc, imin, iminf, iminc, imajh, iminh) = class_fields(indel_keys)

    indel_frame_effect = ""
    if has_indel
        d = length(imaj) - length(ir)
        indel_frame_effect = d % 3 != 0 ? "frameshift" :
                             (d > 0 ? "inframe_insertion" : "inframe_deletion")
    end

    write(snp_fh, join([
        string(location), seq_id, reference_strain, string(is_coding), variant_type,
        string(distinct_strain_count), string(called_strain_count),
        string(no_call_strain_count), call_rate, string(total_weight),
        ref_allele_frequency, string(het_strain_count),
        sr, smaj, smajf, smajc, smin, sminf, sminc, smajh, sminh,
        ir, imaj, imajf, imajc, imin, iminf, iminc, imajh, iminh,
        indel_frame_effect
    ], "\t"), "\n")
end
```

- [ ] **Step 5: Run to verify it passes.** Run the julia command. Expected: PASS for both new testsets; no regressions.

- [ ] **Step 6: Commit.**
```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "$(printf 'variationFeature: per-class SNP/indel schema keyed by (ref,base)\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 3: Rewrite `write_allele_file` — `(ref,base)` keying + per-row `matches_reference`

**Files:**
- Modify: `bin/processSequenceVariations.jl` — `write_allele_file` (~1567–1631)
- Test: `testing/t/handleVariantRecord.jl`

- [ ] **Step 1: Write the failing test.** Append to `testing/t/handleVariantRecord.jl`:

```julia
@testset "write_allele_file keeps deletion and reference as distinct rows" begin
    vars = [
        mkvar(strain="S1", reference="A",   base="G"),
        mkvar(strain="S2", reference="ACA", base="A"),
        mkvar(strain="REF", reference="A",  base="A", matches_reference=1),
    ]
    buf = IOBuffer()
    sid = Dict("S1"=>1, "S2"=>2, "REF"=>3)
    write_allele_file(buf, vars, "LmjF.01", 13850, sid)
    rows = [split(l, '\t') for l in filter(!isempty, split(String(take!(buf)), "\n"))]
    # allele col = 3, matches_reference = 9, genomic_hgvs = 10
    del  = [r for r in rows if r[3] == "A" && r[9] == "0"]
    refr = [r for r in rows if r[3] == "A" && r[9] == "1"]
    @test length(del)  == 1                 # deletion A, matches_reference 0
    @test length(refr) == 1                 # reference A, matches_reference 1
    @test occursin("del", del[1][10])       # deletion has a del g.HGVS
    @test refr[1][10] == "."                # reference row g.HGVS is "."
end
```

- [ ] **Step 2: Run to verify it fails.** Run the julia command. Expected: FAIL — old code collapses the two "A"s into one row (only one row with allele "A").

- [ ] **Step 3: Replace `write_allele_file`.** Replace the entire function (~1567–1631) with:

```julia
function write_allele_file(
    allele_fh::IO,
    variations::Vector{Variation},
    seq_id::String,
    location::Int,
    sample_id_map::Dict{String,Int}
)
    (stats, total_weight) = aggregate_locus_alleles(variations)
    for ((ref, allele), st) in stats
        strain_ids = Set{Int}()
        for strain in st.strains
            sid = get(sample_id_map, strain, 0)
            sid > 0 ? push!(strain_ids, sid) :
                @warn "strain not found in sample_id_map, omitted from strain_ids" strain
        end
        ids_str     = "{" * join(sort(collect(strain_ids)), ",") * "}"
        matches_ref = allele == ref ? 1 : 0
        ghgvs       = matches_ref == 1 ? "." : genomic_hgvs(seq_id, location, ref, allele)
        write(allele_fh, join([
            string(location),
            seq_id,
            allele,
            string(length(st.strains)),
            @sprintf("%.4f", st.weight / total_weight),
            @sprintf("%.2f", st.cov_sum / st.entry_count),
            @sprintf("%.2f", st.pct_sum / st.entry_count),
            ids_str,
            string(matches_ref),
            ghgvs
        ], "\t"), "\n")
    end
end
```

- [ ] **Step 4: Run to verify it passes.** Run the julia command. Expected: PASS; no regressions (the Task-5-of-prior-branch allele.dat sum test still holds).

- [ ] **Step 5: Commit.**
```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "$(printf 'allele.dat: key by (ref,base); per-row matches_reference (QA collapse fix)\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 4: Remove now-dead HGVS helpers

**Files:**
- Modify: `bin/processSequenceVariations.jl` — remove `build_allele_refs` (~1865) and `allele_genomic_hgvs` (~1883) if unreferenced
- Test: `testing/t/handleVariantRecord.jl` (remove any tests of the removed helpers)

- [ ] **Step 1: Confirm they are now unused.**
```bash
grep -n "build_allele_refs\|allele_genomic_hgvs" bin/ testing/ -r
```
Expected: references only in their own definitions and possibly their own unit tests. If any *other* production caller remains, STOP and report — do not remove.

- [ ] **Step 2: Remove the two functions** (`build_allele_refs` and `allele_genomic_hgvs`, with their docstrings) from `bin/processSequenceVariations.jl`. Remove any testset in `testing/t/handleVariantRecord.jl` that calls them.

- [ ] **Step 3: Run the full julia suite.** Run the julia command. Expected: PASS, exit 0, no `UndefVarError`.

- [ ] **Step 4: Commit.**
```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "$(printf 'cleanup: drop build_allele_refs/allele_genomic_hgvs (superseded by direct genomic_hgvs)\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 5: Update e2e assertions to the new schema

**Files:**
- Modify: `testing/t/test_mergeExperiments_e2e.py`

Note: the e2e reads a completed pipeline run's `.dat` files. The existing run at
`/home/jbrestel/dnaseq_test/merge` was produced by the OLD writers, so these tests
can only pass after the pipeline is re-run with this branch. Update the assertions
now so they are ready; flag the re-run in the report.

- [ ] **Step 1: Update the `variationFeature` column count.** In `test_mergeExperiments_e2e.py`, `test_variation_feature_column_count`, change `!= 23` to `!= 31` and the message to "expected 31".

- [ ] **Step 2: Update column-index-based `variationFeature` assertions.** Several tests index `r[...]` against the OLD 23-column layout (e.g. `test_variation_feature_major_allele_nonempty` uses `r[4]`, `test_variation_feature_is_coding_binary` uses `r[13]`, `test_variation_feature_distinct_strain_count_in_range` uses `r[10]`). Re-map or replace them against the new header (0-indexed): `is_coding`=3, `variant_type`=4, `distinct_strain_count`=5, `called_strain_count`=6, `no_call_strain_count`=7, `call_rate`=8, `total_ploidy_count`=9, `ref_allele_frequency`=10, `het_strain_count`=11, `snp_major_allele`=13, `indel_major_allele`=22. For assertions about "major allele" that no longer exist locus-wide, retarget to a class column that is present, or delete the assertion if it no longer maps to a real column. Make each change minimal and named.

- [ ] **Step 3: Add per-class invariants.** Add:
```python
def test_variation_feature_mixed_locus_families_disjoint(work_dirs):
    """A MIXED-locus row populates both families; a SNV/INDEL row leaves the other empty.
    Header-driven so column positions are not hard-coded."""
    path = os.path.join(work_dirs['processSeqVars'], 'variationFeature.dat')
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            r = line.rstrip('\n').split('\t')
            vt = r[idx['variant_type']]
            snp_present   = bool(r[idx['snp_major_allele']])
            indel_present = bool(r[idx['indel_major_allele']])
            if vt == 'SNV':
                assert snp_present and not indel_present
            elif vt == 'INDEL':
                assert indel_present and not snp_present
            elif vt == 'MIXED':
                assert snp_present and indel_present
```

- [ ] **Step 4: Update `allele.dat` assertions.** `allele.dat` stays 10 columns (Task 3 preserves the layout), so `test_allele_dat_column_count` (`!= 10`) is unchanged and should still hold. Confirm `test_allele_dat_frequencies_sum_to_one_per_position` still groups by `(location, seq_id)` — it does; no change. No edits expected here unless a test asserted uniqueness of `(location, allele)` (none does).

- [ ] **Step 5: Commit.**
```bash
git add testing/t/test_mergeExperiments_e2e.py
git commit -m "$(printf 'test: update mergeExperiments e2e to per-class variationFeature schema\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

- [ ] **Step 6: Report the required pipeline re-run.** In the completion report, state that `nextflow run main.nf -entry mergeExperiments -profile mergeExperiments` must be re-run and then `pytest testing/t/test_mergeExperiments_e2e.py --run-dir <launch-dir>` executed to validate end-to-end (the existing run predates this schema change).

---

## Task 6: Update spec status and memory

**Files:**
- Modify: the spec `Status:` line; memory files (outside repo)

- [ ] **Step 1: Mark the spec implemented.** Set the `Status:` line in
`docs/superpowers/specs/2026-07-06-variationfeature-per-class-schema-design.md` to
"Implemented on branch `variationfeature-per-class` — unit-verified; e2e pending a
fresh pipeline run." Commit.

- [ ] **Step 2: Update memory** (`/home/jbrestel/.claude/projects/.../memory/`): in
`project_matches_reference_overlapping_ref_bug.md` and
`project_snp_deletion_allele_collapse.md`, note both are FIXED by the `(ref,base)`
keying on branch `variationfeature-per-class` (pending pipeline re-run). Update
`MEMORY.md` index lines accordingly.

---

## Self-review notes
- **Spec coverage:** 31-col schema (Task 2 header+writer); locus-wide denominator (Task 1 `total`, ref counted once via `:reference` bucket in Task 2); `(ref,base)` identity (Task 1); reference-once + per-class major/minor (Task 2); allele.dat fix (Task 3); dead-helper cleanup (Task 4); testing incl. mixed-locus + e2e (Tasks 2/3/5); both deferred bugs closed (Tasks 2/3, recorded Task 6). Covered.
- **Type consistency:** `AlleleStat` fields (`weight`, `strains`, `cov_sum`, `pct_sum`, `entry_count`) defined Task 1, consumed Tasks 2–3; `aggregate_locus_alleles` returns `(Dict{Tuple{String,String},AlleleStat}, Int)` used identically in both writers; `classify_allele` returns `:reference/:snp/:indel` used in Task 2.
- **Column order** in the Task 2 header string exactly matches the `join([...])` field order in `write_snp_feature` (12 locus-invariant, 9 snp, 10 indel = 31).
- **Known limitation carried forward:** e2e cannot be run in-session (needs a fresh Docker pipeline run); Task 5 updates assertions and Task 5 Step 6 flags the re-run.
