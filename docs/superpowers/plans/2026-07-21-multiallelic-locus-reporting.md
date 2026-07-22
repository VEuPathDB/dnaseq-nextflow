# Multiallelic Locus Reporting Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix two defects at multiallelic loci in `bin/processSequenceVariations.jl`: the output-VCF biallelic split writes other-alt samples as reference, and `variationFeature.dat` major/minor selection ignores the reference allele.

**Architecture:** Two independent, surgical changes to `bin/processSequenceVariations.jl`. Bug 1 changes one lambda in `remap_sample_for_split` (other-alt index → `.` instead of `0`). Bug 2 changes the pool construction and slot-HGVS logic in `write_snp_feature` so the reference competes for major/minor within each variant class. Both are covered by unit tests in `testing/t/handleVariantRecord.jl`.

**Tech Stack:** Julia 1.10, `Test` stdlib. Tests run in the `veupathdb/dnaseqanalysis:latest` Docker image.

**Spec:** `docs/superpowers/specs/2026-07-21-multiallelic-locus-reporting-design.md`

**Test runner (all steps):**
```bash
docker run --rm --pull always -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest \
  julia testing/t/handleVariantRecord.jl
```
Run this from the repo root. The whole file runs in a few seconds; there is no per-testset runner, so "run the test" means run the whole file and read the summary / failure block.

---

## Task 1: Bug 1 — split remap marks other-alt samples missing

**Files:**
- Modify: `bin/processSequenceVariations.jl` (function `remap_sample_for_split`, ~line 1012; docstring ~line 1004)
- Test: `testing/t/handleVariantRecord.jl` (existing testsets at ~line 1610 and ~line 1622)

The existing testsets currently encode the buggy behavior (other alt → `0`). Per TDD we rewrite them to the desired behavior first, watch them fail against the current code, then fix the code.

- [ ] **Step 1: Rewrite the two remap testsets to the desired (missing) behavior**

Replace the testset at `testing/t/handleVariantRecord.jl:1610-1620` with:

```julia
@testset "remap_sample_for_split marks other-alt slots missing, not reference" begin
    # n_orig_alts=2, target_alt_i=2: slot "1" is a NON-target alt → "." (its call
    # lives in the sibling split record); "." stays "."
    @test remap_sample_for_split("1/.", ["GT"], 2, 2) == "./."
    @test remap_sample_for_split("./1", ["GT"], 2, 2) == "./."
    # target_alt_i=1: slot "1" IS the target alt → 1; "." stays "."
    @test remap_sample_for_split("1/.", ["GT"], 2, 1) == "1/."
    # reference slot stays reference
    @test remap_sample_for_split("0/1", ["GT"], 2, 1) == "0/1"
    @test remap_sample_for_split("0/1", ["GT"], 2, 2) == "0/."
    # full-missing GT is left as-is (guarded earlier)
    @test remap_sample_for_split("./.", ["GT"], 2, 2) == "./."
    # sanity: existing biallelic behavior unaffected (n_orig_alts=1 → unchanged)
    @test remap_sample_for_split("1/2", ["GT"], 1, 1) == "1/2"
end
```

Replace the testset at `testing/t/handleVariantRecord.jl:1622-1630` with:

```julia
@testset "remap_sample_for_split remaps triploid GTs, other-alt → missing" begin
    # Triploid split multiallelic. Target alt → 1, reference (0) stays 0,
    # every other alt index → "." ; '.' slots and both separators preserved.
    @test remap_sample_for_split("./2/.", ["GT"], 3, 2) == "./1/."   # target alt → 1
    @test remap_sample_for_split("./2/.", ["GT"], 3, 1) == "././."   # non-target alt → .
    @test remap_sample_for_split("1/./.", ["GT"], 3, 1) == "1/./."   # target alt → 1
    @test remap_sample_for_split("1/2/3", ["GT"], 3, 2) == "./1/."   # only target survives
    @test remap_sample_for_split("0/2/0", ["GT"], 3, 2) == "0/1/0"   # ref slots stay ref
    @test remap_sample_for_split("././.", ["GT"], 3, 2) == "././."   # full-missing untouched
end
```

- [ ] **Step 2: Run the tests to verify the two rewritten testsets FAIL**

Run: the docker command above.
Expected: FAIL. The failure block shows e.g. `remap_sample_for_split("1/.", ["GT"], 2, 2)` produced `"0/."` but expected `"./."` (current code maps other-alt → `0`).

- [ ] **Step 3: Fix the remap lambda**

In `bin/processSequenceVariations.jl`, in `remap_sample_for_split`, replace these lines (~1022-1027):

```julia
            remap = idx -> idx == 0 ? 0 : (idx == target_alt_i ? 1 : 0)
            # Remap every numeric allele index in place: target alt → 1, ref/other → 0.
            # '.' slots (from --multi-overlaps . on a split multiallelic) and the
            # separators (/ or |) are preserved verbatim, so this is ploidy-agnostic
            # — a triploid "./2/." remaps each index without leaving stale ones.
            result[fi] = replace(gt, r"\d+" => m -> string(remap(parse(Int, m))))
```

with:

```julia
            remap = idx -> idx == 0 ? "0" : (idx == target_alt_i ? "1" : ".")
            # Remap every numeric allele index in place: target alt → 1, ref (0) → 0,
            # any OTHER alt → "." (missing) — the sample's real call for that allele
            # lives in the sibling split record, so it must not be fabricated as
            # reference here (the --multi-overlaps . convention). '.' slots and the
            # separators (/ or |) are preserved verbatim, so this is ploidy-agnostic
            # — a triploid "./2/." remaps each index without leaving stale ones.
            result[fi] = replace(gt, r"\d+" => m -> remap(parse(Int, m)))
```

Also update the docstring bullet at ~line 1008. Replace:

```julia
#   - any other alt index → 0 (treated as ref in this split record)
```

with:

```julia
#   - any other alt index → "." (missing; its call lives in the sibling split record)
```

(The docstring lines are `  - target alt index ... → 1`, `  - ref (0) → 0`, `  - any other alt index → 0 ...` — edit the third line's `0 (treated as ref in this split record)` to `"." (missing; its call lives in the sibling split record)`.)

- [ ] **Step 4: Run the tests to verify they PASS**

Run: the docker command above.
Expected: PASS. All testsets in the file green, no errors or warnings.

- [ ] **Step 5: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "Fix multiallelic split: mark other-alt samples missing, not reference

remap_sample_for_split wrote a sample carrying a different ALT as GT=0
(homozygous reference) in each biallelic split row, inflating the
homozygous-reference counts a genome browser reads per record. Map other
alt indices to '.' (missing) instead, matching the --multi-overlaps .
convention already used at the merge step.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: Bug 2 — reference competes for major/minor in write_snp_feature

**Files:**
- Modify: `bin/processSequenceVariations.jl` (function `write_snp_feature`, ~lines 1603-1684)
- Test: `testing/t/handleVariantRecord.jl` (new testsets; plus updates to existing testsets at ~769, ~791, ~1557)

### 2a: Add new failing tests for reference inclusion

- [ ] **Step 1: Add new testsets**

Add these testsets in `testing/t/handleVariantRecord.jl` immediately after the existing `write_snp_feature variant_type MIXED ...` testset (after ~line 927):

```julia
@testset "write_snp_feature: reference wins minor when it outranks a rarer alt" begin
    # Mirrors Pf3D7_01_v3:481838 (haploid): T=2 strains, C=1 strain, ref G=2 strains
    # (incl. synthetic ref). Ranking over {ref,alts}: major=T(2)? tie with ref(2).
    # Give T weight 3 so ordering is unambiguous: major=T, minor=ref G, C dropped.
    vT1 = Variation(); vT1.strain="s1"; vT1.base="T"; vT1.reference="G"; vT1.ploidy=1
    vT2 = Variation(); vT2.strain="s2"; vT2.base="T"; vT2.reference="G"; vT2.ploidy=1
    vT3 = Variation(); vT3.strain="s3"; vT3.base="T"; vT3.reference="G"; vT3.ploidy=1
    vC  = Variation(); vC.strain="s4";  vC.base="C"; vC.reference="G"; vC.ploidy=1
    vR1 = Variation(); vR1.strain="s5"; vR1.base="G"; vR1.reference="G"; vR1.ploidy=1
    vRef= Variation(); vRef.strain="ref"; vRef.base="G"; vRef.reference="G"; vRef.ploidy=1

    buf = IOBuffer()
    write_snp_feature(buf, [vT1,vT2,vT3,vC,vR1,vRef], 1, "Pf3D7_01_v3", 481838, "ref",
                      ["s1","s2","s3","s4","s5"])
    f = split(chomp(String(take!(buf))), '\t')
    @test f[13] == "G"                              # snp_ref_allele
    @test f[14] == "T"                              # snp_major_allele (weight 3)
    @test f[16] == "3"                              # snp_major_allele_strain_count
    @test f[20] == "Pf3D7_01_v3:g.481838G>T"        # snp_major_genomic_hgvs
    @test f[17] == "G"                              # snp_minor_allele = REFERENCE
    @test f[19] == "2"                              # snp_minor_allele_strain_count (s5 + ref)
    @test f[21] == ""                               # snp_minor_genomic_hgvs blank for ref slot
end

@testset "write_snp_feature: reference is major at a typical ref-majority SNP" begin
    # ref C: 3 strains (incl. synthetic ref) ; alt T: 1 strain. major=ref, minor=alt.
    vR1 = Variation(); vR1.strain="s1"; vR1.base="C"; vR1.reference="C"; vR1.ploidy=1
    vR2 = Variation(); vR2.strain="s2"; vR2.base="C"; vR2.reference="C"; vR2.ploidy=1
    vT  = Variation(); vT.strain="s3";  vT.base="T"; vT.reference="C"; vT.ploidy=1
    vRef= Variation(); vRef.strain="ref"; vRef.base="C"; vRef.reference="C"; vRef.ploidy=1

    buf = IOBuffer()
    write_snp_feature(buf, [vR1,vR2,vT,vRef], 0, "chr1", 50, "ref", ["s1","s2","s3"])
    f = split(chomp(String(take!(buf))), '\t')
    @test f[14] == "C"                 # snp_major_allele = reference
    @test f[16] == "3"                 # major strain count (s1,s2,ref)
    @test f[20] == ""                  # major hgvs blank for ref slot
    @test f[17] == "T"                 # snp_minor_allele = the alt
    @test f[19] == "1"                 # minor strain count
    @test f[21] == "chr1:g.50C>T"      # minor hgvs
end

@testset "write_snp_feature: indel reference competes for indel major/minor" begin
    # ref span ACA: 2 strains (incl synthetic ref); deletion A: 1 strain.
    # indel major = reference span ACA (weight 2), minor = deletion A.
    vR1 = Variation(); vR1.strain="s1"; vR1.base="ACA"; vR1.reference="ACA"; vR1.ploidy=1
    vD  = Variation(); vD.strain="s2";  vD.base="A";   vD.reference="ACA"; vD.ploidy=1
    vRef= Variation(); vRef.strain="ref"; vRef.base="ACA"; vRef.reference="ACA"; vRef.ploidy=1

    buf = IOBuffer()
    write_snp_feature(buf, [vR1,vD,vRef], 0, "chr1", 200, "ref", ["s1","s2"])
    f = split(chomp(String(take!(buf))), '\t')
    @test f[22] == "ACA"               # indel_ref_allele
    @test f[23] == "ACA"               # indel_major_allele = reference span
    @test f[29] == ""                  # indel major hgvs blank for ref slot
    @test f[26] == "A"                 # indel_minor_allele = deletion
    @test occursin("del", f[30])       # indel minor hgvs is a deletion
end

@testset "write_snp_feature: class with no alts stays empty (monoallelic)" begin
    # only the reference allele present -> SNP and indel families both empty
    vR1 = Variation(); vR1.strain="s1"; vR1.base="A"; vR1.reference="A"; vR1.ploidy=1
    vRef= Variation(); vRef.strain="ref"; vRef.base="A"; vRef.reference="A"; vRef.ploidy=1
    buf = IOBuffer()
    write_snp_feature(buf, [vR1,vRef], 0, "chr1", 10, "ref", ["s1"])
    f = split(chomp(String(take!(buf))), '\t')
    @test f[14] == ""                  # snp_major_allele empty (no snp alt)
    @test f[17] == ""                  # snp_minor_allele empty
    @test f[23] == ""                  # indel_major_allele empty
end
```

Column indices (1-based) for reference while reading assertions: 13 `snp_ref_allele`, 14 `snp_major_allele`, 15 `snp_major_allele_frequency`, 16 `snp_major_allele_strain_count`, 17 `snp_minor_allele`, 18 `snp_minor_allele_frequency`, 19 `snp_minor_allele_strain_count`, 20 `snp_major_genomic_hgvs`, 21 `snp_minor_genomic_hgvs`, 22 `indel_ref_allele`, 23 `indel_major_allele`, 26 `indel_minor_allele`, 29 `indel_major_genomic_hgvs`, 30 `indel_minor_genomic_hgvs`.

- [ ] **Step 2: Run the tests to verify the new testsets FAIL**

Run: the docker command above.
Expected: FAIL. Current code ranks alts only, so e.g. "reference wins minor" fails at `f[17] == "G"` (current minor is `C`, the rarer alt), and "reference is major" fails at `f[14] == "C"` (current major is `T`).

### 2b: Implement reference-inclusive ranking

- [ ] **Step 3: Rewrite the pool construction and class_fields in write_snp_feature**

In `bin/processSequenceVariations.jl`, replace the block that builds `ref_weight`, `snp_keys`, `indel_keys` and ranks them (currently ~lines 1614-1631):

```julia
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
    # tie-break: weight desc, then allele string, then ref span — the last key
    # keeps ordering deterministic for same-string alleles from different refs
    rank(ks) = sort(ks; by = k -> (-stats[k].weight, k[2], k[1]))
    snp_keys   = rank(snp_keys)
    indel_keys = rank(indel_keys)
```

with:

```julia
    ref_weight   = 0
    snp_alts     = Tuple{String,String}[]
    indel_alts   = Tuple{String,String}[]
    ref_keys     = Tuple{String,String}[]
    for (key, st) in stats
        c = classify_allele(key[1], key[2])
        if c == :reference
            ref_weight += st.weight
            push!(ref_keys, key)
        elseif c == :snp
            push!(snp_alts, key)
        else
            push!(indel_alts, key)
        end
    end
    # The reference competes for a class's major/minor slot, but only when that
    # class actually has an alt allele. Match the reference key to a class by its
    # ref-string (so an MNP's multi-base reference stays in the SNP class, and an
    # indel's ref span stays in the indel class).
    snp_refstrings   = Set(k[1] for k in snp_alts)
    indel_refstrings = Set(k[1] for k in indel_alts)
    snp_keys   = isempty(snp_alts)   ? Tuple{String,String}[] :
                 vcat(snp_alts,   [k for k in ref_keys if k[1] in snp_refstrings])
    indel_keys = isempty(indel_alts) ? Tuple{String,String}[] :
                 vcat(indel_alts, [k for k in ref_keys if k[1] in indel_refstrings])
    # tie-break: weight desc, then allele string, then ref span — deterministic for
    # same-string alleles from different refs; the reference sorts by its own allele.
    rank(ks) = sort(ks; by = k -> (-stats[k].weight, k[2], k[1]))
    snp_keys   = rank(snp_keys)
    indel_keys = rank(indel_keys)
```

Then, in the same function, update `has_snp`/`has_indel` (currently ~lines 1641-1642) so `variant_type` stays keyed on ALT presence, not pool presence:

```julia
    has_snp   = !isempty(snp_alts)
    has_indel = !isempty(indel_alts)
```

(Delete the old `has_snp = !isempty(snp_keys)` / `has_indel = !isempty(indel_keys)` lines.)

- [ ] **Step 4: Blank the genomic_hgvs for a reference slot in class_fields**

In `write_snp_feature`, the nested `class_fields` closure (currently ~lines 1646-1662) computes `majh` and `mnh`. Replace its body:

```julia
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
```

with (only the two `*h` HGVS lines change — blank when the slot holds the reference, i.e. allele == ref):

```julia
    class_fields = function(keys)
        isempty(keys) && return ("","","","","","","","","")
        ref  = keys[1][1]
        maj  = keys[1][2]
        majf = @sprintf("%.4f", stats[keys[1]].weight / total_weight)
        majc = string(length(stats[keys[1]].strains))
        majh = maj == ref ? "" : genomic_hgvs(seq_id, location, ref, maj)
        if length(keys) > 1
            mn   = keys[2][2]
            mnf  = @sprintf("%.4f", stats[keys[2]].weight / total_weight)
            mnc  = string(length(stats[keys[2]].strains))
            mnh  = mn == keys[2][1] ? "" : genomic_hgvs(seq_id, location, keys[2][1], mn)
        else
            mn = ""; mnf = ""; mnc = ""; mnh = ""
        end
        (ref, maj, majf, majc, mn, mnf, mnc, majh, mnh)
    end
```

- [ ] **Step 5: Run the tests to verify the new testsets PASS**

Run: the docker command above.
Expected: the four new testsets PASS. Existing testsets `write_snp_feature emits 31 columns ...`, `write_snp_feature indel major_genomic_hgvs for a deletion`, and `write_snp_feature emits per-class SNP+indel columns without collapse` will now FAIL (they assert the old alt-only major/minor) — fixed in 2c. All other testsets stay green.

### 2c: Update existing tests that asserted alt-only major/minor

- [ ] **Step 6: Update the three affected existing testsets**

In `testing/t/handleVariantRecord.jl`, testset `write_snp_feature emits 31 columns, no CDS fields` (~769-789): the pool is now `{ref A (weight 2: ref+s2), alt T (weight 1: s1)}`, so major = reference A, minor = T. Replace its assertions block (the `@test fields[13]`/`[14]`/`[20]` lines) with:

```julia
    @test fields[13] == "A"                    # snp_ref_allele
    @test fields[14] == "A"                    # snp_major_allele = reference (weight 2)
    @test fields[20] == ""                     # snp_major_genomic_hgvs blank for ref slot
    @test fields[17] == "T"                    # snp_minor_allele = the alt
    @test fields[21] == "LmjF.01:g.500A>T"     # snp_minor_genomic_hgvs
```

In testset `write_snp_feature indel major_genomic_hgvs for a deletion` (~791-804): pool is `{ref CA (weight 2: ref+s2), del C (weight 1: s1)}`, so indel major = reference CA, minor = deletion C. Replace its assertion block with:

```julia
    @test fields[5]  == "INDEL"                    # variant_type
    @test fields[22] == "CA"                       # indel_ref_allele
    @test fields[23] == "CA"                       # indel_major_allele = reference span (weight 2)
    @test fields[29] == ""                         # indel_major_genomic_hgvs blank for ref slot
    @test fields[26] == "C"                        # indel_minor_allele = deletion
    @test fields[30] == "LmjF.01:g.2532delA"       # indel_minor_genomic_hgvs
```

In testset `write_snp_feature emits per-class SNP+indel columns without collapse` (~1557-1577): SNP pool is `{alt G (weight 1: S1), ref A (weight 1: REF)}` — a weight tie broken by allele string, so `"A" < "G"` ⇒ snp_major = reference A, snp_minor = G. Indel pool is `{del A (weight 1: S2)}` only (no indel reference key present, since the synthetic REF variation has ref-string `A`, not `ACA`), so indel_major = deletion A. Replace the assertion block with:

```julia
    @test length(cols) == 31
    @test cols[1] == "13850"
    @test cols[5] == "MIXED"
    @test cols[13] == "A"                       # snp_ref_allele
    @test cols[14] == "A"                       # snp_major_allele = reference (tie, "A"<"G")
    @test cols[20] == ""                        # snp_major_genomic_hgvs blank for ref slot
    @test cols[17] == "G"                       # snp_minor_allele = the SNP alt
    @test cols[21] == "LmjF.01:g.13850A>G"      # snp_minor_genomic_hgvs
    @test cols[22] == "ACA"                     # indel_ref_allele
    @test cols[23] == "A"                       # indel_major_allele = deletion (no indel ref key present)
    @test occursin("del", cols[29])             # indel_major_genomic_hgvs
    @test cols[31] == "frameshift"              # indel_frame_effect (1-3=-2)
```

- [ ] **Step 7: Run the full test file to verify ALL testsets PASS**

Run: the docker command above.
Expected: PASS. Every testset green, output pristine (no errors/warnings).

- [ ] **Step 8: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "variationFeature: rank the reference allele into major/minor

write_snp_feature ranked major/minor among ALT alleles only, mislabeling
the more-common reference at 98% of SNV loci. The reference now competes
for a class's major/minor slot (per-class, matched by ref-string, only when
that class has an alt). A slot holding the reference emits a blank
genomic_hgvs. Aggregate ref_allele_frequency and all count columns are
unchanged.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: Full unit suite + contract docs + memory

**Files:**
- Verify: all of `testing/t/*.jl`, `testing/t/*.py`
- Modify: `docs/superpowers/specs/2026-07-07-complex-variant-frequency-overcount-design.md` (append a note) — only if it documents major/minor; otherwise skip
- Memory: `project_variant_frequency_contract.md`

- [ ] **Step 1: Run the whole Julia + Python unit suite**

Run:
```bash
docker run --rm --pull always -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest bash -c '
  for t in testing/t/*.jl; do julia "$t"; done
  python3 -m pytest testing/t/ -q'
```
Expected: all Julia testsets pass; Python unit tests pass; the 96 `test_mergeExperiments_e2e.py` tests SKIP (no `--run-dir`). If any non-e2e test fails, fix it before proceeding.

- [ ] **Step 2: Update the frequency-contract memory**

Append to `/home/jbrestel/.claude/projects/-home-jbrestel-workspaces-dataLoad-project-home-dnaseq-nextflow/memory/project_variant_frequency_contract.md` a paragraph recording: major/minor now rank the reference allele into each class (top-2 of {ref + alts}, reference joins only when the class has an alt, ties by allele string); a reference-major/minor slot has blank genomic_hgvs; `snp_ref_allele`/`indel_ref_allele` and `ref_allele_frequency` unchanged; the web app must replicate this ranking. Reference the spec/plan dated 2026-07-21.

- [ ] **Step 3: Commit docs/memory**

```bash
git add -A
git commit -m "Docs: record reference-inclusive major/minor contract change

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

- [ ] **Step 4: Flag e2e validation to John (manual)**

The 96 `test_mergeExperiments_e2e.py` tests assert against a real run's output files and skip without `--run-dir`. Any assertions on `snp_major/minor` columns or per-record genotype counts (homozygous-reference) now expect the new semantics. Full validation requires re-running `mergeExperiments` on this branch and pointing the e2e suite at the fresh run dir:
```bash
python3 -m pytest testing/t/test_mergeExperiments_e2e.py --run-dir /path/to/fresh/run
```
This is a heavier, human-initiated step — report to John rather than running it unprompted. If any e2e assertions on these columns are hardcoded to old values, update them to the new semantics as a follow-up.

---

## Notes for the implementer

- Do NOT touch `aggregate_locus_alleles`, `ref_allele_frequency`, or any count column — the fix is confined to ranking/pool construction and split-GT remapping.
- The reference key's `strains` set includes the synthetic reference strain, so a reference major/minor `strain_count` is one higher than a pure sample count. This is intentional and matches `ref_allele_frequency`.
- Bug 1 and Bug 2 are independent; Task 1 and Task 2 can be done in either order. Task 3 runs last.
