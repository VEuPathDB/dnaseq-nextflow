# Separate SNP/indel handling in mergeExperiments (via `-m both`) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop `bcftools merge` from fusing SNP and indel alleles into mixed multiallelic rows, so the annotator classifies each row unambiguously, indels reach snpeff with CANN frame effects, and `downstream_of_frameshift` is surfaced for indels — while keeping `.dat` output locus-based.

**Architecture:** Change `mergeVcfs` to `bcftools merge --merge both` (SNPs merged with SNPs, indels with indels, never across classes → at most one multiallelic SNP record + one multiallelic indel record per start position). Then in `bin/processSequenceVariations.jl`, process the output/CANN path over *both* class-records at a locus (not a single picked record), fix `./.` cross-record reference synthesis so a strain isn't counted twice at a locus, add a `downstream_frameshift` compound effect for indels, and map the new term in `parseSnpEffAnnotations.py`. snpeff runs once, unchanged.

**Tech Stack:** Nextflow DSL2, bcftools 1.19, Julia 1.10 (`Test` stdlib), Python 3 (`pytest`), Perl `Test2::V0` (integration).

**Spec:** `docs/superpowers/specs/2026-07-06-separate-snp-indel-merge-both-design.md`

**Invariant relied on (empirically verified, bcftools 1.19):** `-m both` yields ≤2 records per start position — one multiallelic SNP record, one multiallelic indel record. Different-length indels + a deletion collapse into one indel record with a **common padded REF** (e.g. `ATG → ATTG,ATTTG,A`); `len_diff` is invariant under that padding.

---

## File structure

| File | Responsibility | Action |
|---|---|---|
| `modules/mergeExperiments.nf` | `mergeVcfs` merge flag | Modify (line 44) |
| `bin/processSequenceVariations.jl` | annotator: variation build, CANN, output VCF | Modify (`build_variations_from_record`, `handle_variant_record!`, `build_cann_string`) |
| `bin/parseSnpEffAnnotations.py` | CANN→snpeff.dat effect mapping | Modify (`CANN_EFFECT_MAP`) |
| `testing/t/handleVariantRecord.jl` | Julia unit/characterization tests | Modify (add tests) |
| `testing/t/test_parseSnpEffAnnotations.py` | Python parser tests | Modify (add test) |

---

## Task 1: Switch merge to `-m both`

**Files:**
- Modify: `modules/mergeExperiments.nf:44`
- Test: `testing/t/mergeBoth.t.sh` (new, standalone shell characterization test)

- [ ] **Step 1: Write the failing test**

Create `testing/t/mergeBoth.t.sh`:

```bash
#!/usr/bin/env bash
# Characterization: bcftools merge -m both keeps SNP and indel alleles on
# SEPARATE records at a shared start position (never a mixed multiallelic row).
set -euo pipefail
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT; cd "$tmp"

cat > hdr.txt <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
EOF
mk() { { cat hdr.txt; printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$1";
         printf 'chr1\t100\t.\t%s\t%s\t.\t.\t.\tGT\t1/1\n' "$2" "$3"; } > "$1.vcf"
       bgzip -f "$1.vcf"; tabix -f -p vcf "$1.vcf.gz"; }
mk SA A G       # SNP
mk SB A AT      # insertion

out="$(bcftools merge -m both SA.vcf.gz SB.vcf.gz 2>/dev/null | grep -vc '^#')"
[ "$out" = "2" ] || { echo "FAIL: expected 2 records, got $out"; exit 1; }

# And prove -m all does NOT (documents the bug we are fixing)
allrows="$(bcftools merge -m all SA.vcf.gz SB.vcf.gz 2>/dev/null | grep -vc '^#')"
[ "$allrows" = "1" ] || { echo "FAIL: expected -m all to fuse into 1 row, got $allrows"; exit 1; }

echo "PASS"
```

- [ ] **Step 2: Run test to verify it passes against bcftools (this test asserts a bcftools property, not our code yet)**

Run: `bash testing/t/mergeBoth.t.sh`
Expected: `PASS` (confirms the flag semantics the fix relies on).

- [ ] **Step 3: Make the change**

In `modules/mergeExperiments.nf`, line 44, change:

```
    bcftools merge --merge all -O z -o merged.vcf.gz *.norm.vcf.gz
```

to:

```
    bcftools merge --merge both -O z -o merged.vcf.gz *.norm.vcf.gz
```

- [ ] **Step 4: Verify the module still parses**

Run: `nextflow run main.nf -entry mergeExperiments -profile mergeExperiments -stub-run`
Expected: workflow completes stub run without a DSL parse error (stubs `touch` outputs).

- [ ] **Step 5: Commit**

```bash
git add modules/mergeExperiments.nf testing/t/mergeBoth.t.sh
git commit -m "mergeVcfs: use bcftools merge -m both to keep SNP/indel rows separate"
```

---

## Task 2 (SPIKE): Characterize `./.` reference synthesis across class-records

Goal: prove, with a test, that today a strain carrying a SNP is *also* given a fabricated reference variation from the sibling indel record's `./.` GT when the position is covered. This locks the behavior we must change before restructuring.

**Files:**
- Test: `testing/t/handleVariantRecord.jl` (add testset)

- [ ] **Step 1: Write the characterization test**

Append to `testing/t/handleVariantRecord.jl`:

```julia
# ---------------------------------------------------------------------------
# -m both: cross-record ./. must not fabricate a duplicate reference variation
# ---------------------------------------------------------------------------

@testset "SPIKE: covered ./. on sibling record synthesizes a reference variation" begin
    # Strain S1 carries the SNP (A>G); on the indel record it is ./. .
    # Coverage says the locus is covered, so build_variations_from_record
    # currently synthesizes a REF variation for S1 from the indel record.
    indel = make_vcf_record(ref="ATG", alts=["A"],
                            format_keys=["GT","DP"], sample_data=["./.:0"])
    chrom_cov = Dict("chr1" => [(1, 1000, 30.0)])   # whole contig covered @30x
    vars = build_variations_from_record(indel, ["S1"], Set{String}(), chrom_cov, 1)

    # CURRENT behavior (the thing we will change in Task 3):
    @test length(vars) == 1
    @test vars[1].matches_reference == 1
    @test vars[1].reference == "ATG"
    @test vars[1].base == "ATG"
end
```

- [ ] **Step 2: Run to confirm it PASSES (documents current behavior)**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: PASS. (This is a characterization test of the current, to-be-changed behavior.)

- [ ] **Step 3: Commit the spike**

```bash
git add testing/t/handleVariantRecord.jl
git commit -m "test: characterize cross-record ./. reference synthesis (spike)"
```

---

## Task 3: Suppress cross-record reference synthesis at a mixed locus

A strain that carries a real variant on one class-record must not also be given a synthesized reference variation from the sibling record's covered `./.`. We make reference-synthesis opt-in via a flag, and have `handle_variant_record!` synthesize reference variations **once per locus** (not once per record).

**Files:**
- Modify: `bin/processSequenceVariations.jl` — `build_variations_from_record` (1104–1140), `handle_variant_record!` (2169–2174)
- Test: `testing/t/handleVariantRecord.jl`

- [ ] **Step 1: Write the failing test**

Replace the SPIKE testset body from Task 2 with the target behavior and add a locus-level test. Append:

```julia
@testset "build_variations_from_record: synthesize_ref=false skips covered ./." begin
    indel = make_vcf_record(ref="ATG", alts=["A"],
                            format_keys=["GT","DP"], sample_data=["./.:0"])
    chrom_cov = Dict("chr1" => [(1, 1000, 30.0)])
    vars = build_variations_from_record(indel, ["S1"], Set{String}(), chrom_cov, 1;
                                        synthesize_ref=false)
    @test isempty(vars)
end

@testset "build_variations_from_record: synthesize_ref=true keeps old behavior" begin
    indel = make_vcf_record(ref="ATG", alts=["A"],
                            format_keys=["GT","DP"], sample_data=["./.:0"])
    chrom_cov = Dict("chr1" => [(1, 1000, 30.0)])
    vars = build_variations_from_record(indel, ["S1"], Set{String}(), chrom_cov, 1;
                                        synthesize_ref=true)
    @test length(vars) == 1
    @test vars[1].matches_reference == 1
end
```

Delete the old `SPIKE:` testset (it asserted the pre-change behavior).

- [ ] **Step 2: Run to verify it fails**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: FAIL — `synthesize_ref` is not a keyword argument yet (MethodError / UndefKeywordError).

- [ ] **Step 3: Add the `synthesize_ref` keyword**

In `bin/processSequenceVariations.jl`, change the signature of `build_variations_from_record` (line 1104–1110) to add a keyword:

```julia
function build_variations_from_record(
    record::VCFRecord,
    all_strains::Vector{String},
    undone_strains::Set{String},
    chrom_coverage::Dict{String, Vector{Tuple{Int, Int, Float64}}},
    ploidy::Int=1;
    synthesize_ref::Bool=true
)::Vector{Variation}
```

Then guard the reference-synthesis block. Change the `if covered` block (lines 1122–1138) so it only runs when `synthesize_ref` is true — wrap the existing body:

```julia
        if isempty(gt) || gt == "." || gt == "./." || gt == ".|."
            # No call: synthesize a reference Variation if position is covered
            synthesize_ref || continue
            (covered, dp) = get_coverage(chrom_coverage, strain, record.pos - 1)
            if covered
                # ... unchanged body ...
            end
            continue
        end
```

- [ ] **Step 4: Run to verify the unit tests pass**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: PASS for both new testsets.

- [ ] **Step 5: Make `handle_variant_record!` synthesize reference once per locus**

In `handle_variant_record!` (lines 2169–2174), replace the combined-variation build so records contribute only *real* variations, and reference synthesis happens once against the SNP-preferred record. Change:

```julia
    # Build variations from all records at this position (SNPs + indels combined)
    variations = Variation[]
    for r in records
        append!(variations, build_variations_from_record(r, all_strains, ctx.undone_strains, chrom_coverage, ctx.ploidy))
    end
    isempty(variations) && return false
```

to:

```julia
    # Build variations from all records at this position (SNPs + indels combined).
    # Reference synthesis for covered no-calls happens ONCE per locus (against the
    # SNP-preferred record) so a strain carrying a variant on one class-record is
    # not also given a fabricated reference variation from the sibling record.
    ref_record = pick_snp_record(records)
    variations = Variation[]
    for r in records
        append!(variations, build_variations_from_record(
            r, all_strains, ctx.undone_strains, chrom_coverage, ctx.ploidy;
            synthesize_ref = (r === ref_record)))
    end
    isempty(variations) && return false
```

- [ ] **Step 6: Run the full Julia test file**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: PASS (no regressions in existing testsets).

- [ ] **Step 7: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "annotator: synthesize reference calls once per locus, not per class-record"
```

---

## Task 4: Add `downstream_frameshift` compound effect for indels in CANN

An indel that sits downstream of an upstream frameshift should carry both its structural effect and `downstream_frameshift`. `build_cann_string` currently returns the pure-indel structural effect early (line 1952) with no frameshift context.

**Files:**
- Modify: `bin/processSequenceVariations.jl` — `build_cann_string` (1943–1953)
- Test: `testing/t/handleVariantRecord.jl`

- [ ] **Step 1: Write the failing test**

Append to `testing/t/handleVariantRecord.jl`:

```julia
@testset "build_cann_string: indel downstream of frameshift → compound effect" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=42,
                          pos_in_codon_val=2)
    v = Variation()
    v.downstream_of_frameshift = 1
    # deletion ATG>A : len_diff = -2 → frameshift structurally
    s = build_cann_string("ATG", "A", v, ann)
    @test occursin("frameshift&downstream_frameshift", s)
end

@testset "build_cann_string: indel NOT downstream of frameshift → structural only" begin
    ann = make_annotation(is_coding=1)
    v = Variation()
    v.downstream_of_frameshift = 0
    s = build_cann_string("ATG", "A", v, ann)   # frameshift
    @test occursin("|frameshift|", s)
    @test !occursin("downstream_frameshift", s)
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: FAIL — current output is `k0|.|.|frameshift|...`, missing `&downstream_frameshift`.

- [ ] **Step 3: Implement**

In `build_cann_string`, replace the `is_pure_indel` block (lines 1943–1953):

```julia
    if is_pure_indel
        len_diff = alt_len - ref_len
        structural = if len_diff % 3 != 0
            "frameshift"
        elseif len_diff > 0
            "inframe_insertion"
        else
            "inframe_deletion"
        end
        effect = v.downstream_of_frameshift == 1 ?
            "$(structural)&downstream_frameshift" : structural
        return "k0|.|.|$(effect)|$(tid)|$(pos_in_cds)|$(pic)|.|."
    end
```

- [ ] **Step 4: Run to verify it passes**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: PASS for both new testsets, no regressions.

- [ ] **Step 5: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "CANN: mark indels downstream of a frameshift with compound effect"
```

---

## Task 5: Emit both class-records in the output/CANN path

Today CANN accumulation (line 2238) and the output write loop (2265–2287) run off a single `record = pick_snp_record(records)`, so indel alts never reach `output.vcf`. Restructure so both the SNP record and the indel record contribute CANN entries and output rows, each reading GT from its own sample columns.

**Files:**
- Modify: `bin/processSequenceVariations.jl` — `handle_variant_record!` (2176, 2205–2244, 2262–2287)
- Test: `testing/t/handleVariantRecord.jl`

- [ ] **Step 1: Write the failing test**

Append to `testing/t/handleVariantRecord.jl` a helper that captures written VCF rows and asserts an indel alt is emitted. Add:

```julia
@testset "handle_variant_record!: indel record produces an output.vcf row" begin
    # SNP record: S1 has A>G. Indel record: S2 has ATG>A. Two strains.
    snp   = make_vcf_record(ref="A",   alts=["G"], format_keys=["GT","DP"],
                            sample_data=["1/1:30", "./.:0"])
    indel = make_vcf_record(ref="ATG", alts=["A"], format_keys=["GT","DP"],
                            sample_data=["./.:0", "1/1:30"])
    chrom_cov = Dict("chr1" => [(1, 1000, 30.0)])

    ctx = make_test_context(reference_strain="REF", all_strains=["S1","S2"],
                            ploidy=1)   # see helper below
    io  = IOBuffer()
    writers = make_test_writers(vcf_io=io)   # see helper below

    handle_variant_record!([snp, indel], Tuple{String,Int}[], ctx, writers,
        TranscriptSequenceCache(Dict{String, Dict{String,String}}()),
        ["S1","S2"], chrom_cov)

    out = String(take!(io))
    @test occursin("\tATG\tA\t", out)   # indel row emitted
    @test occursin("\tA\tG\t", out)     # SNP row emitted
end
```

> **Note for implementer:** `make_test_context` and `make_test_writers` are test scaffolding. Build them at the top of the test file to construct a minimal `ProcessingContext` and `OutputWriters` writing to `IOBuffer`s, using an intergenic (`is_coding=0`) position so no transcript DB is needed. If a non-coding position makes `handle_variant_record!` return before writing, use `make_annotation(is_coding=1)` via a stubbed `annotate_position_all` — check the existing `open_output_writers`/`ProcessingContext` constructors and mirror them. If scaffolding proves heavy, split this task: first extract the output-writing section into a testable helper `write_locus_outputs!(...)` and test that directly.

- [ ] **Step 2: Run to verify it fails**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: FAIL — only the SNP row (`A G`) is emitted; the indel row (`ATG A`) is absent.

- [ ] **Step 3: Restructure the record handling**

3a. Keep per-record variation lists so CANN/GT reads come from the correct class-record. After the variation-build block from Task 3, add a parallel structure. Replace the single-record annotation usage.

Change line 2176 region — instead of a single `record = pick_snp_record(records)`, build:

```julia
    # Per-record variation lists: each record's variations read GT from that
    # record's own sample columns. Reference synthesis already happened once
    # (Task 3) against ref_record and lives in that record's list.
    per_record = Tuple{VCFRecord, Vector{Variation}}[]
    for r in records
        rv = build_variations_from_record(
            r, all_strains, ctx.undone_strains, chrom_coverage, ctx.ploidy;
            synthesize_ref = (r === ref_record))
        push!(per_record, (r, rv))
    end
    variations = reduce(vcat, (rv for (_, rv) in per_record); init=Variation[])
    isempty(variations) && return false
```

Remove the now-duplicated variation-build loop introduced in Task 3 Step 5 (this supersedes it — `variations` is derived from `per_record`).

3b. In the annotation loop (2205–2244), change the CANN collection call to iterate every record with its own variations. Replace lines 2238–2243:

```julia
        for (rec, recvars) in per_record
            for (alt, smap) in collect_cann_entries_for_annotation(recvars, annotation, rec, ctx.reference_strain, all_strains, strain_idx_map)
                for (strain, entries) in smap
                    strain_map = get!(alt_strain_entries, alt, Dict{String, Vector{String}}())
                    append!(get!(strain_map, strain, String[]), entries)
                end
            end
        end
```

3c. Replace the output write loop (2262–2287) to iterate every record. Replace:

```julia
    (ref_keys, ref_cann_entries) = build_ref_cann_entries(annotations)
    (alt_cann_entries, alt_strain_to_ca) = assign_cann_keys(alt_strain_entries, all_strains)

    for (rec, recvars) in per_record
        modified_sample_data = fill_missing_coverage_gt(rec, all_strains, chrom_coverage, ctx.ploidy)
        strain_to_dfs = Dict{String,String}(v.strain => string(v.downstream_of_frameshift) for v in recvars)
        dfs_values = [get(strain_to_dfs, s, ".") for s in all_strains]
        n_orig_alts = length(rec.alts)
        for (alt_i, alt) in enumerate(rec.alts)
            if alt == "*"
                @warn "Unexpected * allele in VCF output at $(seq_id):$(location) — should have been removed by mergeVariantsByLocation.py"
                continue
            end
            haskey(alt_cann_entries, alt) || continue
            coding_alt_entries = filter(!=((".")), alt_cann_entries[alt])
            all_entries = [ref_cann_entries..., coding_alt_entries...]
            full_cann   = isempty(all_entries) ? "." : join(all_entries, ',')
            strain_to_alt_keys = get(alt_strain_to_ca, alt, Dict{String, Vector{String}}())
            ca_values = build_ca_values(rec.format_keys, modified_sample_data, rec.alts, alt,
                                        ref_keys, all_strains, strain_to_alt_keys)
            write_vcf_entry(writers.vcf_fh, seq_id, location, rec.ref, alt,
                            full_cann, rec.info, rec.format_keys,
                            modified_sample_data, ca_values, dfs_values,
                            n_orig_alts, alt_i)
        end
    end
```

> **Important:** `alt_strain_entries` is keyed by alt-allele string. SNP alts (`G`,`C`) and indel alts (`A`,`ATTG`) are disjoint strings, so keying by alt across both records does not collide. Verify no two records at a locus share an identical alt string (they cannot: same string with same REF would be the same allele and bcftools would have merged them into one record).

- [ ] **Step 4: Run to verify it passes**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: PASS — both SNP and indel rows emitted; existing testsets green.

- [ ] **Step 5: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "annotator: emit both SNP and indel class-records to output.vcf and CANN"
```

---

## Task 6: Map `downstream_frameshift` in the snpeff.dat parser

The parser must render CANN `downstream_frameshift` (and the compound `frameshift&downstream_frameshift`) into a defined effect/impact. Compound `&` splitting already exists.

**Files:**
- Modify: `bin/parseSnpEffAnnotations.py` — `CANN_EFFECT_MAP` (line 47+)
- Test: `testing/t/test_parseSnpEffAnnotations.py`

- [ ] **Step 1: Write the failing test**

Append to `testing/t/test_parseSnpEffAnnotations.py`:

```python
def test_cann_downstream_frameshift_maps_to_defined_effect():
    from parseSnpEffAnnotations import CANN_EFFECT_MAP
    assert "downstream_frameshift" in CANN_EFFECT_MAP
    effect, impact = CANN_EFFECT_MAP["downstream_frameshift"]
    assert effect == "downstream_of_frameshift_variant"
    assert impact == "MODERATE"
```

- [ ] **Step 2: Run to verify it fails**

Run: `python3 -m pytest testing/t/test_parseSnpEffAnnotations.py::test_cann_downstream_frameshift_maps_to_defined_effect -v`
Expected: FAIL — KeyError / assertion (term not in map).

- [ ] **Step 3: Implement**

In `bin/parseSnpEffAnnotations.py`, add to `CANN_EFFECT_MAP`:

```python
    "downstream_frameshift": ("downstream_of_frameshift_variant", "MODERATE"),
```

- [ ] **Step 4: Run to verify it passes**

Run: `python3 -m pytest testing/t/test_parseSnpEffAnnotations.py -v`
Expected: PASS (new test + no regressions).

- [ ] **Step 5: Commit**

```bash
git add bin/parseSnpEffAnnotations.py testing/t/test_parseSnpEffAnnotations.py
git commit -m "parser: map CANN downstream_frameshift to a defined snpeff.dat effect"
```

---

## Task 7: Integration — run the pipeline test profile and refresh expectations

**Files:**
- Possibly modify: golden expectation files under `testing/` (whatever the test profile diffs against)

- [ ] **Step 1: Run the full Julia + Python unit suites**

Run:
```bash
julia testing/t/handleVariantRecord.jl
python3 -m pytest testing/t/test_parseSnpEffAnnotations.py -v
bash testing/t/mergeBoth.t.sh
```
Expected: all PASS.

- [ ] **Step 2: Run the Nextflow test profile**

Run: `nextflow run main.nf -entry runTests -profile tests`
Expected: Perl `Test2::V0` suite runs. Some `.dat`/VCF golden comparisons for the merge path may now differ because (a) indel rows now appear in `output.vcf`, (b) no duplicate reference variations at mixed loci, (c) indel CANN effects present.

- [ ] **Step 3: Inspect any diffs before accepting**

For each failing golden comparison, confirm the diff matches the intended change:
- indel alleles now present in `output.vcf` / snpeff.dat (source `product_call` and `snpeff`);
- a strain carrying a SNP no longer also emits a reference allele for the sibling indel record;
- `downstream_frameshift` appears only where an indel is downstream of a frameshift.
Do NOT blindly overwrite goldens — verify each change is expected, then update the fixture.

- [ ] **Step 4: Commit refreshed expectations (only after verification)**

```bash
git add testing/
git commit -m "test: refresh merge-path goldens for separated SNP/indel handling"
```

- [ ] **Step 5: Update the spec status**

Mark the spec `Status:` line as `Implemented` and commit:

```bash
git add docs/superpowers/specs/2026-07-06-separate-snp-indel-merge-both-design.md
git commit -m "docs: mark separate-snp-indel spec implemented"
```

---

## Self-review notes

- **Spec coverage:** `-m both` (Task 1); emit indels to output.vcf + one snpeff run (Task 5, snpeff process unchanged); indel CANN frame effects (already in `build_cann_string`, exercised via Task 5); `downstream_frameshift` for indels (Task 4 + parser Task 6); `./.` no-fabrication risk (Tasks 2–3); ≤2-rows invariant + padded REF (relied on in Task 5 alt-key note, fixtures use padded `ATG`). All covered.
- **Type consistency:** `synthesize_ref` keyword used identically in Tasks 3 and 5; `per_record::Vector{Tuple{VCFRecord,Vector{Variation}}}` introduced in Task 5 Step 3a and consumed in 3b/3c; `build_cann_string`/`collect_cann_entries_for_annotation`/`build_ca_values`/`fill_missing_coverage_gt` signatures used as defined in the source.
- **Risk carried into execution:** Task 5 Step 1 scaffolding (`make_test_context`/`make_test_writers`) depends on the real `ProcessingContext`/`OutputWriters` constructors; the step gives a fallback (extract `write_locus_outputs!` and test it directly) if scaffolding is heavy.
- **Out of scope (unchanged):** deferred `matches_reference` overlapping-ref bug; per-experiment `processSingleExperiment` VCF processing.
