# Phase 2: Genomic / Coding / Protein HGVS Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add HGVS `g.` (genomic), `c.` (coding), and `p.` (protein) strings to the three tables whose grain each notation belongs to — `g.`→`allele.dat`, `c.`→`snpeff.dat`, `p.`→`transcript_product.dat`.

**Architecture:** Each notation lives where its coordinate system lives. `g.` is per-allele and transcript-free (a new `genomic_hgvs` reducer handling sub/del/ins/delins from the allele's own `(ref, alt)`); `c.` is per-allele×transcript (pulled through the SnpEff-annotation parser from the `ANN` and `CANN` fields it already reads); `p.` is per-codon (a `protein_hgvs` helper factored out of the existing `substitution_hgvs`). Substitutions are single-nucleotide only for `c.`/`p.` (Phase 1 scope); `g.` covers all allele shapes via the reducer.

**Tech Stack:** Julia 1.10 (`bin/processSequenceVariations.jl`), Julia `Test` (`testing/t/handleVariantRecord.jl`), Python 3 (`bin/parseSnpEffAnnotations.py`, `testing/t/test_parseSnpEffAnnotations.py`). Julia runs in the `jbrestel/dnaseq:latest` container: `docker run --rm -v "$PWD":/work -w /work --entrypoint /opt/julia/bin/julia jbrestel/dnaseq:latest testing/t/handleVariantRecord.jl`.

**Scope / representation notes:**
- `g.` indels use **left-aligned** VCF coordinates (NOT 3'-shifted). Consistent internal ID; not guaranteed to match an external HGVS for repeat-region indels. `dup` is not emitted (plain `ins` is used) — valid HGVS, less idiomatic. This is a deliberate, accepted scope choice.
- `g.` strings are stored **fully qualified** (`{seq_id}:g....`) so record page and results summary read one value, never re-concatenate.
- Deleted bases ARE included in `del`/`delins` (`g.2532delA`, `g.101_102delinsCCTG`) to match SnpEff's `c.` style that appears alongside.
- `c.`: `snpeff`-source rows carry SnpEff's `ANN` HGVS.c verbatim (covers coding + intron + UTR + indels in SnpEff's notation); `product_call`-source rows carry our CANN `hgvs_c` (CDS substitutions only, `.` otherwise). They agree for coding SNVs (verified 284/284 in Phase 1).
- `p.`: consensus-aware, per codon row; reference-codon rows render `p.<AA><pos>=`; frameshift-suppressed rows are already skipped by `write_transcript_product`.
- **Cache safety:** the cache reads only the first 4 columns of `transcript_product.dat` (`seq_id, location, transcript_id, pos_in_cds`); appending `hgvs_p` at the end is safe. `allele.dat` is not a cache.

---

## File Structure

- `bin/processSequenceVariations.jl`
  - **New**: `genomic_hgvs(seq_id, pos, ref, alt) -> String` — pure reducer (Task 1).
  - **New**: `protein_hgvs(protpos, ref_aa, alt_aa) -> String` — factored out of `substitution_hgvs` (Task 3).
  - **Modify**: `substitution_hgvs` to call `protein_hgvs` (Task 3).
  - **Modify**: `write_allele_file` — track per-allele ref, emit `genomic_hgvs` column (Task 2).
  - **Modify**: `write_transcript_product` — emit `protein_hgvs` column (Task 4).
  - **Modify**: `allele.dat` header (line 1253) and `transcript_product.dat` header (line 1255).
- `bin/parseSnpEffAnnotations.py`
  - **Modify**: `parse_ann_rows` and `parse_cann_rows` to yield `hgvs_c`; `main()` header + output row (Task 5).
- `testing/t/handleVariantRecord.jl` — tests for `genomic_hgvs`, `protein_hgvs`, and the two writers.
- `testing/t/test_parseSnpEffAnnotations.py` — tests for the `hgvs_c` column.

---

### Task 1: `genomic_hgvs` reducer

**Files:**
- Modify: `bin/processSequenceVariations.jl` (add near `substitution_hgvs`, after line 1783)
- Test: `testing/t/handleVariantRecord.jl` (new testset near the `substitution_hgvs` tests)

- [ ] **Step 1: Write the failing tests**

```julia
# ---------------------------------------------------------------------------
# genomic_hgvs — per-allele g. (sub / del / ins / delins), left-aligned
# ---------------------------------------------------------------------------

@testset "genomic_hgvs substitution" begin
    @test genomic_hgvs("LmjF.01", 3745, "G", "C") == "LmjF.01:g.3745G>C"
end
@testset "genomic_hgvs substitution embedded in shared prefix" begin
    @test genomic_hgvs("LmjF.01", 100, "CA", "CG") == "LmjF.01:g.101A>G"
end
@testset "genomic_hgvs pure deletion, single base" begin
    @test genomic_hgvs("LmjF.01", 2531, "CA", "C") == "LmjF.01:g.2532delA"
end
@testset "genomic_hgvs pure deletion, multi base" begin
    @test genomic_hgvs("LmjF.01", 3037, "CGA", "C") == "LmjF.01:g.3038_3039delGA"
end
@testset "genomic_hgvs pure insertion" begin
    @test genomic_hgvs("LmjF.01", 1585, "C", "CGA") == "LmjF.01:g.1585_1586insGA"
end
@testset "genomic_hgvs delins from complex allele" begin
    @test genomic_hgvs("LmjF.01", 100, "ATA", "ACCTG") == "LmjF.01:g.101_102delinsCCTG"
end
@testset "genomic_hgvs delins from equal-length MNV" begin
    @test genomic_hgvs("LmjF.01", 100, "CA", "GT") == "LmjF.01:g.100_101delinsGT"
end
@testset "genomic_hgvs uppercases soft-masked bases" begin
    @test genomic_hgvs("LmjF.01", 3745, "g", "c") == "LmjF.01:g.3745G>C"
end
@testset "genomic_hgvs returns dot for no change" begin
    @test genomic_hgvs("LmjF.01", 100, "A", "A") == "."
end
```

- [ ] **Step 2: Run to verify they fail**

Run: `docker run --rm -v "$PWD":/work -w /work --entrypoint /opt/julia/bin/julia jbrestel/dnaseq:latest testing/t/handleVariantRecord.jl 2>&1 | grep -i "genomic_hgvs\|error"`
Expected: `UndefVarError: genomic_hgvs not defined`.

- [ ] **Step 3: Implement**

Add after `substitution_hgvs` (after line 1783):

```julia
"""
    genomic_hgvs(seq_id, pos, ref, alt) -> String

Builds a fully-qualified genomic HGVS string ("{seq_id}:g....") for one allele
from its own (ref, alt) pair at genomic position `pos`. Reduces the pair by
stripping the common prefix/suffix, then classifies: single-base substitution
(">"), pure insertion ("ins"), pure deletion ("del"), or complex/MNV
("delins"). Indels are left-aligned (not 3'-shifted). Returns "." for no change.
"""
function genomic_hgvs(seq_id::String, pos::Int, ref::String, alt::String)::String
    ref = uppercase(ref)
    alt = uppercase(alt)

    # strip common prefix
    p = 0
    while p < min(length(ref), length(alt)) && ref[p+1] == alt[p+1]
        p += 1
    end
    # strip common suffix (without overlapping the stripped prefix)
    s = 0
    while s < min(length(ref), length(alt)) - p && ref[end-s] == alt[end-s]
        s += 1
    end

    r = ref[p+1 : length(ref)-s]
    a = alt[p+1 : length(alt)-s]
    start = pos + p

    if isempty(r) && isempty(a)
        return "."
    elseif length(r) == 1 && length(a) == 1
        return "$(seq_id):g.$(start)$(r)>$(a)"
    elseif isempty(r)                      # insertion between start-1 and start
        return "$(seq_id):g.$(start-1)_$(start)ins$(a)"
    elseif isempty(a)                      # deletion
        return length(r) == 1 ? "$(seq_id):g.$(start)del$(r)" :
                                "$(seq_id):g.$(start)_$(start+length(r)-1)del$(r)"
    else                                   # delins (complex / MNV)
        return length(r) == 1 ? "$(seq_id):g.$(start)delins$(a)" :
                                "$(seq_id):g.$(start)_$(start+length(r)-1)delins$(a)"
    end
end
```

- [ ] **Step 4: Run to verify pass**

Run the full suite; expected exit 0, all 9 `genomic_hgvs` testsets pass, nothing else breaks.

- [ ] **Step 5: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "Add genomic_hgvs reducer (sub/del/ins/delins) for Phase 2 g."
```

---

### Task 2: emit `g.` in `allele.dat`

**Files:**
- Modify: `bin/processSequenceVariations.jl` — `write_allele_file` (lines 1554-1614) and header (line 1253)
- Test: `testing/t/handleVariantRecord.jl`

Context: `write_allele_file` keys `allele_entries` by allele string. `g.` needs each allele's own ref, so we build a parallel `allele_refs: allele -> Set{ref}`. Reference alleles (never seen as an alt) and alleles with an ambiguous multi-ref collision emit `.`. The `make_variation` test helper already exists (grep `function make_variation`); it sets `reference`, `base`, `alt_allele`, `strain`, `coverage`, `percent`, `ploidy`.

- [ ] **Step 1: Write the failing tests**

```julia
@testset "write_allele_file emits genomic_hgvs column for alt alleles" begin
    # ref C at 700; s1 hom alt T; s2 hom ref C
    v_ref = Variation(); v_ref.strain="ref"; v_ref.base="C"; v_ref.reference="C"; v_ref.ploidy=1; v_ref.coverage="10"; v_ref.percent="100"
    v_s1  = Variation(); v_s1.strain="s1";  v_s1.base="T"; v_s1.reference="C"; v_s1.ploidy=1; v_s1.coverage="12"; v_s1.percent="100"
    v_s2  = Variation(); v_s2.strain="s2";  v_s2.base="C"; v_s2.reference="C"; v_s2.ploidy=1; v_s2.coverage="9";  v_s2.percent="100"

    buf = IOBuffer()
    write_allele_file(buf, [v_ref, v_s1, v_s2], "LmjF.01", 700, Dict("ref"=>5,"s1"=>1,"s2"=>2))
    rows = Dict{String,Vector{SubString{String}}}()
    for ln in filter(!isempty, split(String(take!(buf)), "\n"))
        f = split(ln, "\t"); rows[f[3]] = f          # key by allele (col 3)
    end
    @test length(rows["T"]) == 10                      # 9 original columns + genomic_hgvs
    @test rows["T"][10] == "LmjF.01:g.700C>T"          # alt allele gets g.
    @test rows["C"][10] == "."                          # reference allele -> dot
end

@testset "write_allele_file genomic_hgvs for a deletion allele" begin
    v_ref = Variation(); v_ref.strain="ref"; v_ref.base="CA"; v_ref.reference="CA"; v_ref.ploidy=1; v_ref.coverage="10"; v_ref.percent="100"
    v_s1  = Variation(); v_s1.strain="s1";  v_s1.base="C";  v_s1.reference="CA"; v_s1.ploidy=1; v_s1.coverage="12"; v_s1.percent="100"

    buf = IOBuffer()
    write_allele_file(buf, [v_ref, v_s1], "LmjF.01", 2531, Dict("ref"=>5,"s1"=>1))
    rows = Dict{String,Vector{SubString{String}}}()
    for ln in filter(!isempty, split(String(take!(buf)), "\n"))
        f = split(ln, "\t"); rows[f[3]] = f
    end
    @test rows["C"][10] == "LmjF.01:g.2532delA"
end
```

- [ ] **Step 2: Run to verify they fail**

Expected: `BoundsError`/length mismatch — rows currently have 9 columns, `[10]` is absent.

- [ ] **Step 3: Implement**

In `write_allele_file`, after the `allele_entries` build loop (after line 1575), add the per-allele ref map:

```julia
    # per-allele reference span, for genomic HGVS. Only alt alleles get an entry;
    # reference alleles (and multi-ref collisions) fall through to ".".
    allele_refs = Dict{String, Set{String}}()
    for v in variations
        if !isempty(v.alt_allele)
            push!(get!(allele_refs, v.alt_allele, Set{String}()), v.reference)
        elseif v.base != v.reference
            push!(get!(allele_refs, v.base, Set{String}()), v.reference)
        end
    end
```

Then inside the `for (allele, entries)` loop, before the `write(...)`, compute the string:

```julia
        genomic_hgvs_str = (haskey(allele_refs, allele) && length(allele_refs[allele]) == 1) ?
            genomic_hgvs(seq_id, location, first(allele_refs[allele]), allele) : "."
```

And append it to the written row (after `string(matches_ref)`):

```julia
            string(matches_ref),
            genomic_hgvs_str
```

Update the header (line 1253) to append `\tgenomic_hgvs`:

```julia
    write(allele_fh, "location\tseq_id\tallele\tdistinct_strain_count\tallele_frequency\tavg_coverage\tavg_percent\tstrain_ids\tmatches_reference\tgenomic_hgvs\n")
```

- [ ] **Step 4: Run to verify pass** — full suite exit 0, both new testsets pass.

- [ ] **Step 5: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "Emit genomic_hgvs (g.) column in allele.dat"
```

---

### Task 3: extract `protein_hgvs` from `substitution_hgvs`

**Files:**
- Modify: `bin/processSequenceVariations.jl` — add `protein_hgvs`, refactor `substitution_hgvs` (lines 1753-1783)
- Test: `testing/t/handleVariantRecord.jl`

- [ ] **Step 1: Write the failing tests**

```julia
@testset "protein_hgvs missense / synonymous / start-loss / stop-loss / unknown" begin
    @test protein_hgvs(320, "D", "H") == "p.Asp320His"
    @test protein_hgvs(134, "T", "T") == "p.Thr134="
    @test protein_hgvs(1,   "M", "T") == "p.Met1?"
    @test protein_hgvs(34,  "Q", "*") == "p.Gln34Ter"
    @test protein_hgvs(327, "*", "Q") == "."          # stop-loss out of scope
    @test protein_hgvs(10,  "M", "X") == "."          # unknown amino acid
end
```

- [ ] **Step 2: Run to verify they fail** — `UndefVarError: protein_hgvs not defined`.

- [ ] **Step 3: Implement**

Add `protein_hgvs` immediately before `substitution_hgvs` (before line 1745):

```julia
"""
    protein_hgvs(protpos, ref_aa, alt_aa) -> String

Builds an HGVS.p string for a single-residue change at protein position
`protpos`. Returns "." when it cannot form one (unknown amino acid, or stop-loss
which needs extension notation — out of scope).
"""
function protein_hgvs(protpos::Int, ref_aa::String, alt_aa::String)::String
    ref3 = get(AA_THREE_LETTER, ref_aa, "")
    alt3 = get(AA_THREE_LETTER, alt_aa, "")
    (isempty(ref3) || isempty(alt3)) && return "."
    if ref_aa == alt_aa
        return "p.$(ref3)$(protpos)="
    elseif ref_aa == "*"
        return "."
    elseif protpos == 1 && ref_aa == "M"
        return "p.Met1?"
    else
        return "p.$(ref3)$(protpos)$(alt3)"
    end
end
```

Then in `substitution_hgvs`, replace the whole `hgvs_p = "."` … block (lines 1766-1780) with:

```julia
    protpos = div(pos_in_cds - 1, 3) + 1
    hgvs_p  = protein_hgvs(protpos, ref_aa, alt_aa)
```

- [ ] **Step 4: Run to verify pass** — full suite exit 0. The 6 new `protein_hgvs` tests pass AND all 8 existing `substitution_hgvs` tests still pass (behavior unchanged by the extraction).

- [ ] **Step 5: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "Extract protein_hgvs helper from substitution_hgvs (DRY)"
```

---

### Task 4: emit `p.` in `transcript_product.dat`

**Files:**
- Modify: `bin/processSequenceVariations.jl` — `write_transcript_product` (lines 1616-1668) and header (line 1255)
- Test: `testing/t/handleVariantRecord.jl`

- [ ] **Step 1: Write the failing tests**

```julia
@testset "write_transcript_product emits hgvs_p per codon row" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=958,
                          pos_in_codon_val=1, ref_codon="GAC", ref_product="D")
    v1 = make_variation(strain="s1", codon="CAC", product=["H"], downstream_of_frameshift=0)
    buf = IOBuffer()
    write_transcript_product(buf, [v1], ann, "LmjF.01", 3745, Dict{String,Int}())
    fields = split(filter(!isempty, split(String(take!(buf)), "\n"))[1], "\t")
    @test length(fields) == 13          # 12 original + hgvs_p
    @test fields[13] == "p.Asp320His"    # protpos = div(958-1,3)+1 = 320
end

@testset "write_transcript_product hgvs_p is synonymous form for reference codon" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                          pos_in_codon_val=1, ref_codon="ATG", ref_product="M")
    v1 = make_variation(strain="s1", codon="ATG", product=["M"], downstream_of_frameshift=0)
    buf = IOBuffer()
    write_transcript_product(buf, [v1], ann, "chr1", 100, Dict{String,Int}())
    fields = split(filter(!isempty, split(String(take!(buf)), "\n"))[1], "\t")
    @test fields[13] == "p.Met4="        # protpos = div(10-1,3)+1 = 4
end
```

- [ ] **Step 2: Run to verify they fail** — length is 12, `[13]` absent.

- [ ] **Step 3: Implement**

In `write_transcript_product`, inside the `for ec in seen_codons` loop, after `matches_ref_product = ...` (line 1652), add:

```julia
        hgvs_p = protein_hgvs(pos_in_protein, annotation.ref_product, product)
```

Append it to the written row (after `dfs_str`):

```julia
            dfs_str,
            hgvs_p
```

Update the header (line 1255) to append `\thgvs_p`:

```julia
    write(tp_fh, "#seq_id\tlocation\ttranscript_id\tpos_in_cds\tpos_in_protein\tcodon\tpos_in_codon\tcount\tproduct\tmatches_ref_codon\tmatches_ref_product\tdownstream_of_frameshift_strain_ids\thgvs_p\n")
```

- [ ] **Step 4: Run to verify pass** — full suite exit 0; both new tests pass; existing `write_transcript_product` tests still pass (they assert `length(fields) == 12`, so **update those to `== 13`** — there are assertions at the existing "12 columns" testset and the DFS testset; change `12`→`13` where a full-row length is asserted, and `fields[12]`→still `fields[12]` for `dfs_str` since it stays column 12). Verify by rereading those testsets before editing.

- [ ] **Step 5: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "Emit hgvs_p (p.) column in transcript_product.dat"
```

---

### Task 5: add `hgvs_c` (c.) to `snpeff.dat`

**Files:**
- Modify: `bin/parseSnpEffAnnotations.py` — `parse_ann_rows` (lines 96-112), `parse_cann_rows` (lines 131-144), `main()` (lines 152, 164-178)
- Test: `testing/t/test_parseSnpEffAnnotations.py`

- [ ] **Step 1: Write the failing tests** — inspect the existing test file first (`python3 -m pytest testing/t/test_parseSnpEffAnnotations.py -q` to see structure), then add tests asserting the new `hgvs_c` column. Add:

```python
def test_ann_row_carries_hgvs_c():
    info = "ANN=G|missense_variant|MODERATE|GENE|GENE|transcript|T1:mRNA|protein_coding|1/1|c.958G>C|p.Asp320His|958/999|958/999|320/332||"
    rows = list(parse_ann_rows(info))
    assert rows == [("G", "T1:mRNA", "MODERATE", "missense_variant", "c.958G>C")]

def test_ann_row_hgvs_c_dot_when_absent():
    info = "ANN=G|intergenic_region|MODIFIER|A-B|A-B|intergenic_region|A-B|||n.233C>G||||||"
    rows = list(parse_ann_rows(info))
    # intergenic: feature_type != transcript, HGVS.c field is n. (still captured verbatim if present) or "."
    assert rows and rows[0][1] == ""          # transcript_id empty
    assert rows[0][4] in (".", "n.233C>G")     # accept dot or the n. string per implementation

def test_cann_row_carries_hgvs_c():
    info = "CANN=k0|CAC|H|missense|T1:mRNA|958|1|c.958G>C|p.Asp320His"
    rows = list(parse_cann_rows("C", info))
    assert rows == [("C", "T1:mRNA", "MODERATE", "missense_variant", "c.958G>C")]
```

> Note: adjust the expected `impact`/`effect` mapping in the CANN test to whatever `map_cann_effect("missense", "CAC")` actually returns — run it to confirm before locking the assertion.

- [ ] **Step 2: Run to verify they fail** — the yielded tuples have 4 elements, not 5.

- [ ] **Step 3: Implement**

In `parse_ann_rows`, capture HGVS.c (ANN field index 9) and add it to the yield. After `transcript_id = ...` (line 106) add:

```python
        hgvs_c = parts[9] if len(parts) > 9 and parts[9] else "."
```
and change the yield (line 112) to:
```python
            yield (allele, transcript_id, impact, component, hgvs_c)
```

In `parse_cann_rows`, capture CANN field index 7. After `transcript_id = parts[4]` (line 138) add:
```python
        hgvs_c = parts[7] if len(parts) > 7 and parts[7] else "."
```
and change the yield (line 144) to:
```python
            yield (alt, transcript_id, impact, so_effect, hgvs_c)
```

In `main()`: change the header (line 152) to add `hgvs_c`:
```python
        out.write("location\tseq_id\tallele\ttranscript_id\timpact\teffect\thgvs_c\tsource\n")
```
Change the unpack + dedup + write (lines 170-178):
```python
            for source, (allele, transcript_id, impact, effect, hgvs_c) in sourced_rows:
                key = (location, seq_id, allele, transcript_id, effect, source)
                if key in seen:
                    continue
                seen.add(key)
                out.write(f"{location}\t{seq_id}\t{allele}\t{transcript_id}\t{impact}\t{effect}\t{hgvs_c}\t{source}\n")
```

- [ ] **Step 4: Run to verify pass**

Run: `python3 -m pytest testing/t/test_parseSnpEffAnnotations.py -q`
Expected: all pass (existing 18 + new). If existing tests assert exact output rows, update them for the new `hgvs_c` column (reread and fix expected strings).

- [ ] **Step 5: Commit**

```bash
git add bin/parseSnpEffAnnotations.py testing/t/test_parseSnpEffAnnotations.py
git commit -m "Add hgvs_c (c.) column to snpeff.dat from ANN and CANN"
```

---

### Task 6: End-to-end regeneration + verification

**Files:** none (verification only).

- [ ] **Step 1: Regenerate allele.dat + transcript_product.dat on real data**

```bash
scratch=$(mktemp -d)
wd=/home/jbrestel/dnaseq_test/merge/work/10/06d80118ee3eee9b998325a951416e
for f in merged.vcf.gz cache.tsv undoneStrains.txt codingSequences.db codingIndels.db lmajF_chr1.gtf coverage.tsv; do cp -L "$wd/$f" "$scratch/$f"; done
docker run --rm -v "$scratch":/work -v "$PWD/bin":/opt/bin -w /work \
  --entrypoint /opt/julia/bin/julia jbrestel/dnaseq:latest /opt/bin/processSequenceVariations.jl \
  --vcf_file merged.vcf.gz --cache_file cache.tsv --undone_strains_file undoneStrains.txt \
  --reference_strain lmajFriedlin --transcript_db codingSequences.db --indel_db codingIndels.db \
  --gtf_file lmajF_chr1.gtf --coverage_file coverage.tsv --ploidy 2 --output_vcf output.vcf 2>&1 | tail -3
echo "SCRATCH=$scratch"
```

- [ ] **Step 2: Verify `g.` forms in allele.dat**

```bash
echo "substitution g.:"; awk -F'\t' 'NR>1 && $10 ~ /g\.[0-9]+[ACGT]>[ACGT]$/' "$scratch/allele.dat" | head -3
echo "deletion g.:";     awk -F'\t' 'NR>1 && $10 ~ /del/'    "$scratch/allele.dat" | head -3
echo "insertion g.:";    awk -F'\t' 'NR>1 && $10 ~ /ins/'    "$scratch/allele.dat" | head -3
echo "delins g.:";       awk -F'\t' 'NR>1 && $10 ~ /delins/' "$scratch/allele.dat" | head -3
echo "reference-allele rows must be '.':"; awk -F'\t' 'NR>1 && $9==1 && $10!="."' "$scratch/allele.dat" | head
```
Expected: SNV/del/ins present; `delins` present at complex/MNV loci; every `matches_reference==1` row has `g.` == `.` (no output from the last check).

- [ ] **Step 3: Verify `p.` in transcript_product.dat**

```bash
awk -F'\t' 'NR>1{print $13}' "$scratch/transcript_product.dat" | grep -c '^p\.'    # populated count
awk -F'\t' 'NR>1 && $13 !~ /^(p\.|\.)$/' "$scratch/transcript_product.dat" | head    # malformed (want none)
```
Expected: p. strings populated; no malformed values.

- [ ] **Step 4: Regenerate snpeff.dat via the parser and verify `hgvs_c`**

```bash
python3 bin/parseSnpEffAnnotations.py --vcf /home/jbrestel/dnaseq_test/merge/output/merged.ann.vcf.gz --output "$scratch/snpeff.dat"
head -1 "$scratch/snpeff.dat"
echo "product_call c. == our CANN c.? (spot)"; awk -F'\t' 'NR>1 && $8=="product_call" && $7 ~ /^c\./' "$scratch/snpeff.dat" | head -3
echo "snpeff-source c. present for coding:"; awk -F'\t' 'NR>1 && $8=="snpeff" && $7 ~ /^c\./' "$scratch/snpeff.dat" | head -3
```
(Confirm the parser CLI flags with `python3 bin/parseSnpEffAnnotations.py --help` first; adjust `--vcf/--output` names if needed.)

- [ ] **Step 5: Cross-check `c.` agreement (product_call vs SnpEff), reuse Phase-1 method**

```bash
awk -F'\t' 'NR>1 && $7 ~ /^c\./{print $7}' "$scratch/snpeff.dat" | sort -u > /tmp/vc.txt
zcat /home/jbrestel/dnaseq_test/merge/output/merged.ann.vcf.gz | grep -v "^#" \
  | grep -oE 'c\.[0-9]+[ACGT]>[ACGT]' | sort -u > /tmp/se.txt
echo "our coding-SNV c. not in SnpEff (want empty):"; comm -23 <(grep -E 'c\.[0-9]+[ACGT]>[ACGT]$' /tmp/vc.txt) /tmp/se.txt | head
rm -f /tmp/vc.txt /tmp/se.txt
```
Expected: empty (every substitution `c.` we emit is also SnpEff's).

- [ ] **Step 6: Clean up** — `rm -rf "$scratch"`.

---

## Self-Review

**Spec coverage:** `g.`→Task 1+2; `p.` helper→Task 3, emitted→Task 4; `c.`→Task 5; validation→Task 6. Header updates in Tasks 2/4/5. Cache safety asserted (Task 4 note). Per-allele ref + collision handling (Task 2).

**Placeholder scan:** none — all steps carry code/commands. Task 5's CANN test asserts a mapped effect; flagged to confirm `map_cann_effect` output before locking.

**Type consistency:** `genomic_hgvs(seq_id::String, pos::Int, ref::String, alt::String)::String` — called in Task 2 with `(seq_id, location, first(allele_refs[allele]), allele)` (String, Int, String, String). `protein_hgvs(protpos::Int, ref_aa::String, alt_aa::String)::String` — called from `substitution_hgvs` (Task 3) and `write_transcript_product` (Task 4) with matching types. Parser yields become 5-tuples consistently across both `parse_ann_rows` and `parse_cann_rows`, unpacked as 5 in `main()`.

**Column counts:** allele.dat 9→10; transcript_product.dat 12→13; snpeff.dat 7→8 (hgvs_c inserted before source).

## Open decisions (defaults chosen)
1. `del`/`delins` include the deleted/replaced bases (`g.2532delA`) to match SnpEff's `c.` style. Drop the trailing bases if you prefer the terse HGVS-2016 form (`g.2532del`).
2. `dup` is not detected; tandem duplications render as `ins`. Valid HGVS; add `dup` detection later if desired.
3. `g.` stored fully qualified with `{seq_id}:` prefix. If you'd rather store the bare `g....` and prepend `seq_id` at display time, drop the prefix in `genomic_hgvs`.
