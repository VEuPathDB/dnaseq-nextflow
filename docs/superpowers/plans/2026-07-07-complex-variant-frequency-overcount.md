# Complex/Multiallelic Variant Frequency Over-Counting — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix allele-frequency over-counting at loci where one strain's diploid genotype is spread across multiple VCF records — complex-variant decomposition (denominator inflation) and `1/2` multiallelic splitting (fabricated reference alleles).

**Architecture:** Two coupled changes. (1) In `mergeVcfs`, tell `bcftools norm -m -any` to mark the complementary allele of a split multiallelic as missing (`.`) instead of reference (`0`). (2) In `processSequenceVariations.jl`, count allele frequencies per **non-missing chromosome slot** and use a **distinct-strain ploidy denominator**, honoring `.` slots as missing. The two changes must land together (see spec Coupling table).

**Tech Stack:** Julia 1.10 (`bin/processSequenceVariations.jl`), bcftools 1.19 (`veupathdb/dnaseqanalysis:1.1.0`), Nextflow DSL2, Perl/Julia test suites.

**Spec:** `docs/superpowers/specs/2026-07-07-complex-variant-frequency-overcount-design.md`

---

## File Map

| File | Change |
|---|---|
| `bin/processSequenceVariations.jl` | Add `allele_slots` field to `Variation` + `chromosome_alleles` helper; harden `gt_to_base` / `nonref_alt_alleles` / `compute_percent` for `.` slots; populate `allele_slots` in `build_variations_from_record`; rewrite `aggregate_locus_alleles` denominator + numerator; fix `het_strain_count` in `write_snp_feature` |
| `modules/mergeExperiments.nf` | Add `--multi-overlaps .` to `bcftools norm -m -any` in `mergeVcfs` |
| `testing/t/handleVariantRecord.jl` | New unit tests for `chromosome_alleles`, `.`-slot parsing, complex/`1/2`-split aggregation, distinct-strain denominator, het_strain_count |

**Design note (backward compatibility):** `chromosome_alleles(v)` returns `v.allele_slots` when populated, else derives the chromosome-copy list from the legacy `base`/`alt_allele`/`ploidy` fields. Existing tests construct `Variation`s with legacy fields only, so they keep passing unchanged; only `build_variations_from_record` sets `allele_slots` (for the new `.`-slot / complex cases).

**Run all Julia unit tests with:** `julia testing/t/handleVariantRecord.jl`

---

### Task 1: Add `allele_slots` field and `chromosome_alleles` helper

**Files:**
- Modify: `bin/processSequenceVariations.jl:450-479` (struct + constructor)
- Modify: `bin/processSequenceVariations.jl` (add helper near `aggregate_locus_alleles`, ~line 1460)
- Test: `testing/t/handleVariantRecord.jl`

- [ ] **Step 1: Write the failing test**

Add near the aggregate tests (after line 1316):

```julia
@testset "chromosome_alleles: legacy hom derives from base×ploidy" begin
    v = Variation(); v.reference = "A"; v.base = "G"; v.ploidy = 2
    @test chromosome_alleles(v) == ["G", "G"]
end

@testset "chromosome_alleles: legacy het derives [ref, alt]" begin
    v = Variation(); v.reference = "A"; v.base = "R"; v.alt_allele = "G"; v.ploidy = 2
    @test chromosome_alleles(v) == ["A", "G"]
end

@testset "chromosome_alleles: explicit allele_slots wins over legacy fields" begin
    v = Variation(); v.reference = "T"; v.base = "TA"; v.ploidy = 2
    v.allele_slots = ["TA"]
    @test chromosome_alleles(v) == ["TA"]
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: FAIL — `UndefVarError: chromosome_alleles not defined` (and `allele_slots` not a field).

- [ ] **Step 3: Add the struct field**

In `bin/processSequenceVariations.jl`, add a field at the end of the `Variation` struct (after `ploidy::Int` on line 473):

```julia
    ploidy::Int
    allele_slots::Vector{String}   # resolved allele per non-missing GT slot; empty ⇒ derive from legacy fields
end
```

Update the `Variation()` constructor (line 476-479) to pass the new default as the last positional argument:

```julia
function Variation()
    Variation("", 0, "", "", "", "", "", "", "", "", "",
              0, 0, 0, 0, "", String[], "", "", 0, 0, 0, 1, String[])
end
```

- [ ] **Step 4: Add the helper**

Insert immediately above `function aggregate_locus_alleles` (currently line 1468):

```julia
"""
    chromosome_alleles(v) -> Vector{String}

The resolved allele string carried by each non-missing chromosome copy this
variation represents. Uses `v.allele_slots` when populated (set by
build_variations_from_record); otherwise derives from the legacy fields: a het
(`alt_allele` set) is `[reference, alt_allele]`; a hom/ref call is `base`
repeated `ploidy` times.
"""
function chromosome_alleles(v::Variation)::Vector{String}
    isempty(v.allele_slots) || return v.allele_slots
    isempty(v.alt_allele) ? fill(v.base, v.ploidy) : [v.reference, v.alt_allele]
end
```

- [ ] **Step 5: Run test to verify it passes**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: PASS (new testsets green; all pre-existing testsets still green).

- [ ] **Step 6: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "feat: add Variation.allele_slots + chromosome_alleles helper"
```

---

### Task 2: Honor `.` (missing) slots in GT parsing and populate `allele_slots`

Under `--multi-overlaps .`, a split `1/2` yields `1/.` and `./1`. Today `gt_to_base` and `nonref_alt_alleles` return empty for any half-missing GT, and `compute_percent` crashes on `parse(Int, ".")`; `build_variations_from_record` then skips the record entirely. Fix all three to skip only the `.` slot, and populate `allele_slots`.

**Files:**
- Modify: `bin/processSequenceVariations.jl:1002-1033` (`gt_to_base`)
- Modify: `bin/processSequenceVariations.jl:1042-1061` (`nonref_alt_alleles`)
- Modify: `bin/processSequenceVariations.jl:1070-1095` (`compute_percent`)
- Modify: `bin/processSequenceVariations.jl:1104-1183` (`build_variations_from_record`)
- Test: `testing/t/handleVariantRecord.jl`

- [ ] **Step 1: Write the failing tests**

Add after the tests from Task 1:

```julia
@testset "gt_to_base: half-missing 1/. returns the present alt" begin
    @test gt_to_base("1/.", "T", ["TA"]) == "TA"
    @test gt_to_base("./1", "T", ["TA"]) == "TA"
end

@testset "nonref_alt_alleles: half-missing keeps the present alt" begin
    @test nonref_alt_alleles("1/.", ["TA"]) == ["TA"]
    @test nonref_alt_alleles("./1", ["TA"]) == ["TA"]
end

@testset "compute_percent: half-missing uses present alt AO" begin
    fmt = Dict("AO" => "7", "RO" => "0")
    @test compute_percent(fmt, "1/.") == "100.00"
    @test compute_percent(fmt, "./1") == "100.00"
end

@testset "build_variations_from_record: 1/. yields one alt slot, no ref" begin
    rec = make_vcf_record(pos=8962, ref="T", alts=["TA"],
                          format_keys=["GT","AO","RO"], sample_data=["1/.:7:0"])
    cov = make_coverage("s1", 8961, 200, 12.0)
    vars = build_variations_from_record(rec, ["s1"], Set{String}(), cov, 2)
    @test length(vars) == 1
    @test vars[1].allele_slots == ["TA"]
    @test vars[1].ploidy == 2
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: FAIL — `gt_to_base("1/.", …)` returns `""`; `compute_percent` throws `ArgumentError` parsing `"."`; `build_variations_from_record` returns an empty vector.

- [ ] **Step 3: Fix `gt_to_base`**

Replace the diploid branch of `gt_to_base` (lines 1011-1032) — the `else` block after `if isnothing(sep_idx)` — with:

```julia
    else
        a1_str = gt[1:sep_idx-1]
        a2_str = gt[sep_idx+1:end]

        # Half-missing (e.g. "1/." from a split multiallelic): return the present allele.
        if a1_str == "." || a2_str == "."
            present = a1_str == "." ? a2_str : a1_str
            present == "." && return ""
            idx = parse(Int, present)
            return idx == 0 ? ref : alts[idx]
        end

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
```

- [ ] **Step 4: Fix `nonref_alt_alleles`**

Replace the diploid index extraction (lines 1047-1051) — the `else` branch of the `idxs = if isnothing(sep_idx)` assignment — with a version that skips `.`:

```julia
    idxs = if isnothing(sep_idx)
        [parse(Int, gt)]
    else
        a1 = gt[1:sep_idx-1]; a2 = gt[sep_idx+1:end]
        collect(parse(Int, s) for s in (a1, a2) if s != ".")
    end
```

- [ ] **Step 5: Fix `compute_percent`**

Replace the diploid branch (lines 1079-1084) — the `else` after `if isnothing(sep_idx)` — with:

```julia
    else
        a1s = gt[1:sep_idx-1]; a2s = gt[sep_idx+1:end]
        a1 = a1s == "." ? 0 : parse(Int, a1s)
        a2 = a2s == "." ? 0 : parse(Int, a2s)
        aidx = a1 != 0 ? a1 : a2
        aidx == 0 && return "0.0"
    end
```

- [ ] **Step 6: Add a slot-resolver and populate `allele_slots` in `build_variations_from_record`**

Insert this helper immediately above `gt_to_base` (before line 1002):

```julia
"""
    resolve_gt_slots(gt, ref, alts) -> Vector{String}

Resolved allele string for each non-missing GT slot. A `.` slot is skipped
(its allele lives in a sibling split-multiallelic record); `0` → ref; `n` → alts[n].
"""
function resolve_gt_slots(gt::String, ref::String, alts::Vector{String})::Vector{String}
    result = String[]
    sep_idx = findfirst(c -> c == '/' || c == '|', gt)
    slots = isnothing(sep_idx) ? [gt] : split(gt, r"[/|]")
    for s in slots
        s == "." && continue
        idx = parse(Int, s)
        push!(result, idx == 0 ? ref : alts[idx])
    end
    result
end
```

In `build_variations_from_record`, set `allele_slots` in **both** places a `Variation` is built:

In the synthesized-reference block, after `v.ploidy = ploidy` (line 1140):

```julia
                v.ploidy             = ploidy
                v.allele_slots       = fill(record.ref, ploidy)
```

In the GT block, after `v.ploidy = gt_ploidy` (line 1177):

```julia
        v.ploidy             = gt_ploidy
        v.allele_slots       = resolve_gt_slots(gt, record.ref, record.alts)
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: PASS (new testsets green; all pre-existing testsets still green — legacy callers don't pass `.` GTs and `chromosome_alleles` still derives from legacy fields for `Variation`s that don't set `allele_slots`).

- [ ] **Step 8: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "fix: honor missing (.) GT slots; populate Variation.allele_slots"
```

---

### Task 3: Rewrite `aggregate_locus_alleles` — per-slot numerator, distinct-strain denominator

**Files:**
- Modify: `bin/processSequenceVariations.jl:1468-1492` (`aggregate_locus_alleles`)
- Test: `testing/t/handleVariantRecord.jl`

- [ ] **Step 1: Write the failing tests**

Add after the existing aggregate test (after line 1316):

```julia
@testset "aggregate: complex decomposition counts each strain's ploidy once" begin
    # LmjF.01:85879 shape — 2 diploid strains carry BOTH a SNP (T>C) and a
    # deletion (TGT>T) via two 1/1 records; 2 ref strains + 1 haploid synthetic ref.
    mkslots(strain, ref, base, ploidy) = begin
        v = Variation(); v.strain = strain; v.reference = ref; v.base = base
        v.ploidy = ploidy; v.percent = "100"; v.coverage = "10"
        v.allele_slots = fill(base, ploidy); v
    end
    vars = [
        mkslots("LV39",     "T",   "C", 2), mkslots("LV39",     "TGT", "T", 2),
        mkslots("LV39cl5",  "T",   "C", 2), mkslots("LV39cl5",  "TGT", "T", 2),
        mkslots("Fried",    "T",   "T", 2), mkslots("Seid",     "T",   "T", 2),
        mkslots("synthref", "T",   "T", 1),
    ]
    (stats, total) = aggregate_locus_alleles(vars)
    @test total == 9                         # 2+2+2+2+1, each strain once
    @test stats[("T","C")].weight   == 4     # both chromosomes of 2 strains
    @test stats[("TGT","T")].weight == 4
    @test stats[("T","T")].weight   == 5     # Fried 2 + Seid 2 + synthref 1
end

@testset "aggregate: split 1/2 compound het fabricates no reference" begin
    # LmjF.01:8962 shape — Seidman is 1/2 (TA/TAA), split into 1/. and ./1.
    het(strain, ref, alt) = begin
        v = Variation(); v.strain = strain; v.reference = ref; v.base = alt
        v.ploidy = 2; v.percent = "100"; v.coverage = "12"
        v.allele_slots = [alt]; v          # one non-missing slot, no ref
    end
    refstrain(strain, ploidy) = begin
        v = Variation(); v.strain = strain; v.reference = "T"; v.base = "T"
        v.ploidy = ploidy; v.percent = "100"; v.coverage = "20"
        v.allele_slots = fill("T", ploidy); v
    end
    vars = [
        het("Seid", "T", "TA"), het("Seid", "T", "TAA"),
        refstrain("LV39", 2), refstrain("Fried", 2), refstrain("LV39cl5", 2),
        refstrain("synthref", 1),
    ]
    (stats, total) = aggregate_locus_alleles(vars)
    @test total == 9                          # Seid counted once as ploidy 2
    @test stats[("T","TA")].weight  == 1
    @test stats[("T","TAA")].weight == 1
    @test stats[("T","T")].weight   == 7      # Seid contributes 0 reference
    @test !("Seid" in stats[("T","T")].strains)
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: FAIL — current code returns `total == 13` / `11` (per-record sum) and puts `"Seid"` in the reference strain set with `("T","T").weight == 9`.

- [ ] **Step 3: Rewrite `aggregate_locus_alleles`**

Replace the whole function body (lines 1468-1492) with:

```julia
function aggregate_locus_alleles(variations::Vector{Variation})::Tuple{Dict{Tuple{String,String},AlleleStat}, Int}
    stats = Dict{Tuple{String,String},AlleleStat}()
    strain_ploidy = Dict{String,Int}()   # each strain's chromosome count, once
    add! = function(ref, allele, weight, strain, cov, pct)
        st = get!(stats, (ref, allele), AlleleStat())
        st.weight      += weight
        push!(st.strains, strain)
        st.cov_sum     += cov
        st.pct_sum     += pct
        st.entry_count += 1
    end
    for v in variations
        strain_ploidy[v.strain] = max(get(strain_ploidy, v.strain, 0), v.ploidy)
        cov     = isempty(v.coverage) ? 0.0 : parse(Float64, v.coverage)
        altfrac = isempty(v.percent)  ? 0.0 : parse(Float64, v.percent)
        slots   = chromosome_alleles(v)
        has_alt = any(a -> a != v.reference, slots)
        for a in slots
            pct = a == v.reference ? (has_alt ? 100.0 - altfrac : altfrac) : altfrac
            add!(v.reference, a, 1, v.strain, cov, pct)
        end
    end
    total = sum(values(strain_ploidy); init=0)
    (stats, total)
end
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: PASS. In particular the pre-existing aggregate/write_allele_file testsets stay green: single-record-per-strain loci have identical totals, and per-slot iteration preserves each allele's weight, strain set, and averaged coverage/percent.

- [ ] **Step 5: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "fix: aggregate_locus_alleles per-slot numerator + distinct-strain denominator"
```

---

### Task 4: Count `het_strain_count` by distinct strain, not by record

**Files:**
- Modify: `bin/processSequenceVariations.jl:1531` (inside `write_snp_feature`)
- Test: `testing/t/handleVariantRecord.jl`

- [ ] **Step 1: Write the failing test**

Add after the aggregate tests:

```julia
@testset "write_snp_feature het_strain_count counts distinct strains" begin
    # One strain, two het records (split 1/2). het_strain_count must be 1, not 2.
    v1 = Variation(); v1.strain="Seid"; v1.reference="T"; v1.base="TA";  v1.alt_allele="TA";  v1.ploidy=2; v1.coverage="12"; v1.percent="58"
    v2 = Variation(); v2.strain="Seid"; v2.reference="T"; v2.base="TAA"; v2.alt_allele="TAA"; v2.ploidy=2; v2.coverage="12"; v2.percent="42"
    buf = IOBuffer()
    write_snp_feature(buf, [v1, v2], 0, "LmjF.01", 8962, "synthref", ["Seid"])
    fields = split(strip(String(take!(buf))), "\t")
    @test fields[12] == "1"   # het_strain_count (column 12)
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: FAIL — `het_strain_count` is `2` (counts variation records).

- [ ] **Step 3: Fix the count**

In `write_snp_feature`, change line 1531 from:

```julia
    het_strain_count = count(v -> !isempty(v.alt_allele), variations)
```

to:

```julia
    het_strain_count = length(Set(v.strain for v in variations if !isempty(v.alt_allele)))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia testing/t/handleVariantRecord.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "fix: het_strain_count counts distinct strains, not records"
```

---

### Task 5: Emit missing (not reference) for split multiallelics in `mergeVcfs`

**Files:**
- Modify: `modules/mergeExperiments.nf:41`

- [ ] **Step 1: Characterization test — confirm the default fabricates reference**

Run (host or container bcftools 1.19):

```bash
printf '##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t100\t.\tT\tTA,TAA\t50\t.\t.\tGT:AD\t1/2:0,7,5\n' > /tmp/mo.vcf
echo "default:";              bcftools norm -m -any               /tmp/mo.vcf 2>/dev/null | grep -v '^##' | cut -f2,4,5,10
echo "--multi-overlaps .:";   bcftools norm -m -any --multi-overlaps . /tmp/mo.vcf 2>/dev/null | grep -v '^##' | cut -f2,4,5,10
```

Expected:
- default → `1/0` (TA) and `0/1` (TAA)  ← fabricated reference
- `--multi-overlaps .` → `1/.` (TA) and `./1` (TAA)  ← correct

- [ ] **Step 2: Edit the process script**

In `modules/mergeExperiments.nf`, change line 41 from:

```
        bcftools norm -m -any -Oz -o "\${vcf%.vcf.gz}.norm.vcf.gz"
```

to:

```
        bcftools norm -m -any --multi-overlaps . -Oz -o "\${vcf%.vcf.gz}.norm.vcf.gz"
```

- [ ] **Step 3: Commit**

```bash
git add modules/mergeExperiments.nf
git commit -m "fix: mergeVcfs uses --multi-overlaps . so split 1/2 does not fabricate reference"
```

---

### Task 6: End-to-end verification on the test dataset

**Files:** none (verification only)

- [ ] **Step 1: Re-run mergeExperiments on the test dataset**

```bash
cd /home/jbrestel/dnaseq_test/merge
nextflow -C nextflow.config run /home/jbrestel/workspaces/dataLoad/project_home/dnaseq-nextflow/main.nf \
  -entry mergeExperiments -profile mergeExperiments -resume
```

Expected: run completes `OK`.

- [ ] **Step 2: Assert the complex locus (85879) is corrected**

```bash
awk -F'\t' '$1==85879' output/variationFeature.dat
```

Expected: `total_ploidy_count` column = `9`, `ref_allele_frequency` = `0.5556`, snp `C` freq = `0.4444`, indel `del` freq = `0.4444`.

- [ ] **Step 3: Assert the compound-het locus (8962) is corrected**

```bash
awk -F'\t' '$1==8962' output/variationFeature.dat
awk -F'\t' '$1==8962' output/allele.dat
```

Expected in `variationFeature.dat`: `total_ploidy_count` = `9`, `ref_allele_frequency` = `0.7778`, `het_strain_count` = `1`, `TA`/`TAA` freq = `0.1111`.
Expected in `allele.dat`: reference (`T`) row strain set is `{1,2,4,5}` (Seidman = strain 3 **absent**); `TA` and `TAA` rows each freq `0.1111`.

- [ ] **Step 4: Spot-check downstream tolerance of `1/.` genotypes**

Confirm `snpEff` and consensus/`CANN` steps did not error on `1/.` / `./1` records:

```bash
nextflow log $(nextflow log | tail -n1 | cut -f3) -f process,exit,workdir | grep -iE "snpEff|processSeqVars|Coding"
zcat output/merged.ann.vcf.gz | awk -F'\t' '$2==8962{print $10,$11,$12,$13}'
```

Expected: all listed processes `exit=0`; the `8962` records show Seidman with `1/.` / `./1` (or the pipeline's post-merge equivalent), no crash.

- [ ] **Step 5: Full unit-test sweep**

```bash
julia testing/t/handleVariantRecord.jl
nextflow run /home/jbrestel/workspaces/dataLoad/project_home/dnaseq-nextflow/main.nf -entry runTests -profile tests
```

Expected: all Julia testsets pass; Perl `Test2::V0` suite passes.

- [ ] **Step 6: Commit any fixture updates**

If the mergeExperiments e2e fixture (`testing/`/spec expected files) encodes 85879/8962 values, update them to the corrected numbers above and commit:

```bash
git add -A
git commit -m "test: update e2e expected values for corrected complex/1-2 frequencies"
```

---

## Self-Review

- **Spec coverage:** Change 1 (merge `--multi-overlaps .`) → Task 5. Change 2a (`.`-slot GT parsing, no record drop) → Task 2. Change 2b (per-slot numerator + distinct-strain denominator) → Tasks 1+3. Change 2c (het_strain_count distinct) → Task 4. Frequency contract → enforced by Tasks 3 (denominator/numerator) and 2 (missing-slot handling); verified in Task 6. VCF-is-correct (no VCF change) → respected (no FreeBayes/atomization edits).
- **Placeholder scan:** none — every code step shows full code; every run step shows the command and expected output.
- **Type consistency:** `chromosome_alleles(v)::Vector{String}` defined in Task 1 and consumed in Task 3; `allele_slots::Vector{String}` field added in Task 1, populated in Task 2, read in Task 1's helper; `resolve_gt_slots` defined and called in Task 2. `het_strain_count` column index (12) matches the header written at line 1255.
