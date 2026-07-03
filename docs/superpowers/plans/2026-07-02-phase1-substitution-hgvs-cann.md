# Phase 1: Substitution HGVS in CANN Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add HGVS `c.` (coding-DNA) and `p.` (protein) nomenclature strings to the CANN annotation for single-base coding substitutions only — no indels, no normalization.

**Architecture:** A new pure helper `substitution_hgvs` builds the two strings from data already computed in `build_cann_string` (CDS position, strand-oriented ref/alt codons, ref/alt amino acids). The CANN entry format gains two trailing fields (`|hgvs_c|hgvs_p`); every entry emits them, but only coding-substitution alt (`k`) entries populate them — all other entries emit `.`. Because the trailing fields are appended and the Python parser reads CANN by fixed low indices with a `len < 5` guard, no functional parser change is required.

**Tech Stack:** Julia 1.10 (`bin/processSequenceVariations.jl`), Julia `Test` (`testing/t/handleVariantRecord.jl`), run in the `jbrestel/dnaseq:latest` container. Python parser docstring/header touch-ups only (`bin/parseSnpEffAnnotations.py`).

**Scope boundary (explicit):** Substitutions only. Indels (`inframe_insertion`/`inframe_deletion`/`frameshift`), complex variants (indel+SNP at anchor), downstream-of-frameshift entries, and ambiguous codons emit `hgvs_c=. hgvs_p=.`. Indel HGVS is out of scope because (a) the pipeline does not compute an indel protein product today, and (b) indel calls are deliberately not 3'-normalized. See `docs/.../` discussion; that is Phase 2+. Stop-loss substitutions (ref amino acid `*` → any residue) also emit `hgvs_p=.` because correct HGVS requires protein-extension notation (`p.Ter327GlnextTer?`); that is deferred to Phase 2+.

**Why the codon is the source of the c. bases:** `ref_codon` and the alt `codon` (`v.codon`) are already built in transcript orientation (strand-correct). The base this variant changed sits at `pos_in_codon` (`pic`). Reading `ref_codon[pic]` and `codon[pic]` therefore yields the strand-correct ref>alt bases with no separate reverse-complement step. `pos_in_cds` is likewise a strand-correct CDS coordinate.

---

## File Structure

- `bin/processSequenceVariations.jl`
  - **New**: `AA_THREE_LETTER` constant (1→3 letter amino-acid map incl. `*`→`Ter`).
  - **New**: `substitution_hgvs(pos_in_cds, ref_codon, alt_codon, pic, ref_aa, alt_aa) -> (hgvs_c, hgvs_p)` — pure, unit-tested in isolation.
  - **Modify**: `build_cann_string` (~line 1760) — call the helper on the substitution branch; append `|hgvs_c|hgvs_p` to every return.
  - **Modify**: `build_ref_cann_entry` (~line 1747) — append `|.|.`.
  - **Modify**: CANN INFO `##INFO=<ID=CANN...>` header string (~line 930) — document the two new fields.
- `testing/t/handleVariantRecord.jl`
  - **New** tests for `substitution_hgvs` and for the extended `build_cann_string`/`build_ref_cann_entry` field count.
- `bin/parseSnpEffAnnotations.py`
  - **Modify**: docstring format comment (~line 121) only. No logic change.

---

### Task 1: `substitution_hgvs` helper + amino-acid table

**Files:**
- Modify: `bin/processSequenceVariations.jl` (add above `# CANN annotation` section header, ~line 1737)
- Test: `testing/t/handleVariantRecord.jl` (add a new testset block after the existing `build_ref_cann_entries` tests, near line 100)

- [ ] **Step 1: Write the failing tests**

Add to `testing/t/handleVariantRecord.jl`:

```julia
# ---------------------------------------------------------------------------
# substitution_hgvs — Phase 1 HGVS for coding substitutions
# ---------------------------------------------------------------------------

@testset "substitution_hgvs missense: c. and p. from strand-oriented codons" begin
    # pos_in_cds=958, codon pos 1 changed G->C, ref Asp(D) -> alt His(H)
    (c, p) = substitution_hgvs(958, "GAC", "CAC", 1, "D", "H")
    @test c == "c.958G>C"
    @test p == "p.Asp320His"           # protpos = div(958-1,3)+1 = 320
end

@testset "substitution_hgvs synonymous uses '=' protein form" begin
    (c, p) = substitution_hgvs(402, "ACT", "ACC", 3, "T", "T")
    @test c == "c.402T>C"
    @test p == "p.Thr134="
end

@testset "substitution_hgvs nonsense uses Ter" begin
    # ref Gln(Q) -> alt stop(*)
    (c, p) = substitution_hgvs(100, "CAA", "TAA", 1, "Q", "*")
    @test c == "c.100C>T"
    @test p == "p.Gln34Ter"            # protpos = div(100-1,3)+1 = 34
end

@testset "substitution_hgvs start-loss renders p.Met1?" begin
    (c, p) = substitution_hgvs(1, "ATG", "ACG", 2, "M", "T")
    @test c == "c.1T>C"
    @test p == "p.Met1?"
end

@testset "substitution_hgvs returns dots for ambiguous codon" begin
    (c, p) = substitution_hgvs(10, "ATG", "ANG", 2, "M", "X")
    @test c == "."
    @test p == "."
end

@testset "substitution_hgvs returns dot p. for unknown amino acid" begin
    (c, p) = substitution_hgvs(10, "ATG", "ACG", 2, "M", "X")
    @test c == "c.10T>C"
    @test p == "."
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run:
```bash
docker run --rm -v "$PWD":/work -w /work --entrypoint /opt/julia/bin/julia \
  jbrestel/dnaseq:latest testing/t/handleVariantRecord.jl 2>&1 | grep -i "substitution_hgvs\|error"
```
Expected: `UndefVarError: substitution_hgvs not defined` (feature missing).

- [ ] **Step 3: Write minimal implementation**

Add to `bin/processSequenceVariations.jl` immediately above the `# CANN annotation` comment banner (~line 1737):

```julia
const AA_THREE_LETTER = Dict(
    "A"=>"Ala","R"=>"Arg","N"=>"Asn","D"=>"Asp","C"=>"Cys",
    "Q"=>"Gln","E"=>"Glu","G"=>"Gly","H"=>"His","I"=>"Ile",
    "L"=>"Leu","K"=>"Lys","M"=>"Met","F"=>"Phe","P"=>"Pro",
    "S"=>"Ser","T"=>"Thr","W"=>"Trp","Y"=>"Tyr","V"=>"Val",
    "*"=>"Ter",
)

"""
    substitution_hgvs(pos_in_cds, ref_codon, alt_codon, pic, ref_aa, alt_aa) -> (String, String)

Builds (HGVS.c, HGVS.p) for a single-base coding substitution. Bases are read
from the strand-oriented codons at position `pic`, so no separate strand
handling is required. Returns "." for either string it cannot form (ambiguous
or non-triplet codon, unknown amino acid, or no base change at `pic`).
"""
function substitution_hgvs(pos_in_cds::Int, ref_codon::String, alt_codon::String,
                           pic::Int, ref_aa::String, alt_aa::String)::Tuple{String,String}
    hgvs_c = "."
    if length(ref_codon) == 3 && length(alt_codon) == 3 &&
       !occursin(r"[^ACGTacgt]", ref_codon) && !occursin(r"[^ACGTacgt]", alt_codon) &&
       1 <= pic <= 3
        rb = uppercase(ref_codon[pic])
        ab = uppercase(alt_codon[pic])
        if rb != ab
            hgvs_c = "c.$(pos_in_cds)$(rb)>$(ab)"
        end
    end

    hgvs_p = "."
    ref3 = get(AA_THREE_LETTER, ref_aa, "")
    alt3 = get(AA_THREE_LETTER, alt_aa, "")
    if !isempty(ref3) && !isempty(alt3)
        protpos = div(pos_in_cds - 1, 3) + 1
        hgvs_p = if ref_aa == alt_aa
            "p.$(ref3)$(protpos)="
        elseif ref_aa == "*"
            "."                       # stop-loss needs extension notation; out of Phase 1 scope
        elseif protpos == 1 && ref_aa == "M"
            "p.Met1?"
        else
            "p.$(ref3)$(protpos)$(alt3)"
        end
    end

    (hgvs_c, hgvs_p)
end
```

- [ ] **Step 4: Run tests to verify they pass**

Run:
```bash
docker run --rm -v "$PWD":/work -w /work --entrypoint /opt/julia/bin/julia \
  jbrestel/dnaseq:latest testing/t/handleVariantRecord.jl 2>&1 | grep -i "substitution_hgvs"
```
Expected: each `substitution_hgvs ...` testset reports all-Pass, and overall exit 0.

- [ ] **Step 5: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "Add substitution_hgvs helper for Phase 1 coding-substitution HGVS"
```

---

### Task 2: Wire HGVS fields into CANN entries

**Files:**
- Modify: `bin/processSequenceVariations.jl` — `build_cann_string` (~line 1760) and `build_ref_cann_entry` (~line 1747)
- Test: `testing/t/handleVariantRecord.jl` (add testset near the existing CANN tests, ~line 100)

- [ ] **Step 1: Write the failing tests**

Add to `testing/t/handleVariantRecord.jl`:

```julia
@testset "build_cann_string appends hgvs_c and hgvs_p for a coding substitution" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=958,
                          pos_in_codon_val=1, ref_codon="GAC", ref_product="D")
    v   = make_variation(strain="s1", codon="CAC", product=["H"])
    entry = build_cann_string("G", "C", v, ann)     # SNP: ref len 1, alt len 1
    parts = split(entry, "|")
    @test length(parts) == 9
    @test parts[8] == "c.958G>C"
    @test parts[9] == "p.Asp320His"
end

@testset "build_cann_string emits dot hgvs for a pure indel" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                          pos_in_codon_val=1, ref_codon="ATG", ref_product="M")
    v   = make_variation(strain="s1", codon=".", product=String[])
    entry = build_cann_string("AT", "A", v, ann)    # deletion
    parts = split(entry, "|")
    @test length(parts) == 9
    @test parts[8] == "."
    @test parts[9] == "."
end

@testset "build_ref_cann_entry appends dot hgvs fields" begin
    ann = make_annotation(is_coding=1, transcript_id="T1", pos_in_cds=10,
                          pos_in_codon_val=1, ref_codon="ATG", ref_product="M")
    entry = build_ref_cann_entry("r0", ann)
    parts = split(entry, "|")
    @test length(parts) == 9
    @test parts[8] == "."
    @test parts[9] == "."
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run:
```bash
docker run --rm -v "$PWD":/work -w /work --entrypoint /opt/julia/bin/julia \
  jbrestel/dnaseq:latest testing/t/handleVariantRecord.jl 2>&1 | grep -i "appends\|emits dot\|error"
```
Expected: FAIL — `length(parts) == 9` fails (currently 7 fields).

- [ ] **Step 3: Write minimal implementation**

In `bin/processSequenceVariations.jl`, edit `build_ref_cann_entry`. Change its final line:

```julia
    "$(key)|$(codon)|$(aa)|reference|$(annotation.transcript_id)|$(annotation.pos_in_cds)|$(annotation.pos_in_codon_val)"
```
to:
```julia
    "$(key)|$(codon)|$(aa)|reference|$(annotation.transcript_id)|$(annotation.pos_in_cds)|$(annotation.pos_in_codon_val)|.|."
```

In `build_cann_string`, append `|.|.` to the pure-indel return:
```julia
        return "k0|.|.|$(structural)|$(tid)|$(pos_in_cds)|$(pic)"
```
→
```julia
        return "k0|.|.|$(structural)|$(tid)|$(pos_in_cds)|$(pic)|.|."
```

Append `|.|.` to the downstream-frameshift return:
```julia
        return "k0|.|.|downstream_frameshift|$(tid)|$(pos_in_cds)|$(pic)"
```
→
```julia
        return "k0|.|.|downstream_frameshift|$(tid)|$(pos_in_cds)|$(pic)|.|."
```

Append `|.|.` to the ambiguous-codon return:
```julia
        return "k0|$(codon)|.|.|$(tid)|$(pos_in_cds)|$(pic)"
```
→
```julia
        return "k0|$(codon)|.|.|$(tid)|$(pos_in_cds)|$(pic)|.|."
```

Replace the substitution (`!is_indel`) return:
```julia
    if !is_indel
        return "k0|$(codon)|$(product_str)|$(aa_effect)|$(tid)|$(pos_in_cds)|$(pic)"
    else
```
→
```julia
    if !is_indel
        alt_aa = length(unique_prods) == 1 ? unique_prods[1] : ""
        (hgvs_c, hgvs_p) = substitution_hgvs(pos_in_cds, annotation.ref_codon, codon, pic,
                                             annotation.ref_product, alt_aa)
        return "k0|$(codon)|$(product_str)|$(aa_effect)|$(tid)|$(pos_in_cds)|$(pic)|$(hgvs_c)|$(hgvs_p)"
    else
```

Append `|.|.` to the complex return (the `else` branch after the substitution branch):
```julia
        return "k0|$(codon)|$(product_str)|$(aa_effect)&$(structural)|$(tid)|$(pos_in_cds)|$(pic)"
```
→
```julia
        return "k0|$(codon)|$(product_str)|$(aa_effect)&$(structural)|$(tid)|$(pos_in_cds)|$(pic)|.|."
```

> Note: the non-coding `annotation.is_coding != 1 && return "."` guards in both functions stay unchanged — a whole-entry `.` has no `|` fields and the Python parser treats it as absent.

- [ ] **Step 4: Run tests to verify they pass**

Run the full Julia suite:
```bash
docker run --rm -v "$PWD":/work -w /work --entrypoint /opt/julia/bin/julia \
  jbrestel/dnaseq:latest testing/t/handleVariantRecord.jl > /tmp/jl.out 2>&1; echo "exit=$?"; grep -ciE "fail|error" /tmp/jl.out
```
Expected: `exit=0` and `0`. (Existing CANN tests use `startswith`/`contains`, so appended fields do not break them.)

- [ ] **Step 5: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "Emit hgvs_c/hgvs_p in CANN entries for coding substitutions"
```

---

### Task 3: Update CANN header + Python parser docstring

**Files:**
- Modify: `bin/processSequenceVariations.jl` (~line 930, the `##INFO=<ID=CANN...>` write)
- Modify: `bin/parseSnpEffAnnotations.py` (~line 121 docstring comment)

- [ ] **Step 1: Update the CANN VCF header description**

In `bin/processSequenceVariations.jl`, the CANN `##INFO` line currently ends:
```
... Format per entry: key|codon|aa|effect|transcript_id|pos_in_cds|pos_in_codon. Compound effects use '&' separator (e.g. missense&frameshift).">
```
Change the format description to:
```
... Format per entry: key|codon|aa|effect|transcript_id|pos_in_cds|pos_in_codon|hgvs_c|hgvs_p. hgvs_c/hgvs_p are populated for coding substitutions only; '.' otherwise. Compound effects use '&' separator (e.g. missense&frameshift).">
```

- [ ] **Step 2: Update the Python parser docstring**

In `bin/parseSnpEffAnnotations.py`, `parse_cann_rows` docstring line:
```python
    CANN entry format: key|codon|aa|effect|transcript_id|pos_in_cds|pos_in_codon
```
→
```python
    CANN entry format: key|codon|aa|effect|transcript_id|pos_in_cds|pos_in_codon|hgvs_c|hgvs_p
```

- [ ] **Step 3: Verify the Python parser still passes (no logic change)**

Run:
```bash
python3 -m pytest testing/t/test_parseSnpEffAnnotations.py -v
```
Expected: PASS (parser reads `parts[0..4]` with a `len(parts) < 5` guard; extra trailing fields are ignored).

- [ ] **Step 4: Commit**

```bash
git add bin/processSequenceVariations.jl bin/parseSnpEffAnnotations.py
git commit -m "Document hgvs_c/hgvs_p CANN fields in header and parser"
```

---

### Task 4: End-to-end regeneration + cross-validation against SnpEff

**Files:** none (verification only). Uses the last merge run's `processSeqVars` inputs.

- [ ] **Step 1: Regenerate the CANN-annotated VCF on real data**

```bash
scratch=$(mktemp -d)
wd=/home/jbrestel/dnaseq_test/merge/work/10/06d80118ee3eee9b998325a951416e
for f in merged.vcf.gz cache.tsv undoneStrains.txt codingSequences.db codingIndels.db lmajF_chr1.gtf coverage.tsv; do cp -L "$wd/$f" "$scratch/$f"; done
docker run --rm -v "$scratch":/work -v "$PWD/bin":/opt/bin -w /work \
  --entrypoint /opt/julia/bin/julia jbrestel/dnaseq:latest /opt/bin/processSequenceVariations.jl \
  --vcf_file merged.vcf.gz --cache_file cache.tsv --undone_strains_file undoneStrains.txt \
  --reference_strain lmajFriedlin --transcript_db codingSequences.db --indel_db codingIndels.db \
  --gtf_file lmajF_chr1.gtf --coverage_file coverage.tsv --ploidy 2 --output_vcf output.vcf
```
Expected: `[ Info: Processed 1000 variant positions`, `output.vcf` written.

- [ ] **Step 2: Confirm coding-substitution k-entries carry HGVS**

```bash
grep -v "^#" "$scratch/output.vcf" | grep -oE 'CANN=[^;	]*' | tr ',' '\n' \
  | grep -E '^k[0-9]+\|[ACGT]{3}\|[A-Z*]\|(missense|synonymous|nonsense)\|' \
  | grep -E 'c\.[0-9]+[ACGT]>[ACGT]\|p\.' | head
```
Expected: substitution k-entries end with populated `...|c.<pos><ref>><alt>|p.<...>`, e.g. `c.958G>C|p.Asp320His`.

- [ ] **Step 3: Cross-validate our c. against SnpEff's ANN c. for coding SNVs**

For each coding SNV, our `hgvs_c` must equal SnpEff's `HGVS.c` (both are reference-based at the DNA level). Compare against the SnpEff-annotated VCF already published (`/home/jbrestel/dnaseq_test/merge/output/merged.ann.vcf.gz`, `ANN` field column 10, `c.` token):

```bash
# our c. strings (from regenerated CANN k-entries), sorted unique
grep -v "^#" "$scratch/output.vcf" | grep -oE 'c\.[0-9]+[ACGT]>[ACGT]' | sort -u > /tmp/ours.txt
# SnpEff c. strings for SNVs
zcat /home/jbrestel/dnaseq_test/merge/output/merged.ann.vcf.gz | grep -v "^#" \
  | grep -oE 'c\.[0-9]+[ACGT]>[ACGT]' | sort -u > /tmp/snpeff.txt
echo "in ours only:"; comm -23 /tmp/ours.txt /tmp/snpeff.txt | head
echo "in snpeff only:"; comm -13 /tmp/ours.txt /tmp/snpeff.txt | head
```
Expected: no `c.` string appears in "ours only". (SnpEff may have extra `c.` for transcripts/positions we treat differently, so "snpeff only" may be non-empty — that is acceptable; investigate only if "ours only" is non-empty, which would indicate a wrong coordinate or base.)

- [ ] **Step 4: Sanity-check p. divergence is confined to multi-SNP codons**

Our `p.` should match SnpEff's `p.` **except** where our consensus codon carries a co-occurring SNP (the intended codon-aware difference). Spot-check any mismatch:
```bash
grep -v "^#" "$scratch/output.vcf" | grep -oE 'c\.[0-9]+[ACGT]>[ACGT]\|p\.[A-Za-z0-9=?]+' | sort -u | head -20
```
Expected: forms look like `c.402T>C|p.Thr134=` (synonymous), `c.958G>C|p.Asp320His` (missense), `c.100C>T|p.Gln34Ter` (nonsense). No crashes, no malformed strings.

- [ ] **Step 5: Clean up scratch**

```bash
rm -rf "$scratch" /tmp/ours.txt /tmp/snpeff.txt /tmp/jl.out
```

---

## Self-Review

**Spec coverage:**
- c./p. for substitutions → Task 1 (`substitution_hgvs`) + Task 2 (wiring).
- Strand correctness → handled implicitly by reading strand-oriented `ref_codon`/`codon` at `pic` (documented in header + Task 1 docstring); validated in Task 4 Step 3.
- Synonymous (`=`), nonsense (`Ter`), missense, start-loss (`Met1?`) → Task 1 tests.
- No normalization / no indel HGVS → enforced by only populating the `!is_indel` branch; asserted in Task 2 ("emits dot hgvs for a pure indel").
- CANN format change documented → Task 3.
- Parser safety → Task 3 Step 3 (regression run); rationale in Architecture.
- Real-data validation + SnpEff cross-check → Task 4.

**Placeholder scan:** none — all steps carry exact code/commands.

**Type consistency:** `substitution_hgvs(pos_in_cds::Int, ref_codon::String, alt_codon::String, pic::Int, ref_aa::String, alt_aa::String)::Tuple{String,String}` is called in Task 2 with `(pos_in_cds, annotation.ref_codon, codon, pic, annotation.ref_product, alt_aa)` — all types match (`codon`/`ref_product` are `String`, `pos_in_cds`/`pic` are `Int`). CANN entry field count is 9 consistently across all returns and all Task 2 assertions.

## Open decisions (confirm before Task 1, low-stakes defaults chosen)
1. **Synonymous form:** plan uses `p.Thr134=` (HGVS-preferred) rather than `p.Thr134Thr`. Change the helper's synonymous branch if the record page prefers the repeated-AA form.
2. **Field placement:** HGVS appended as fields 8–9 of the CANN entry. Alternative would be a separate `HGVSC=`/`HGVSP=` INFO tag; rejected because you asked for it *in* CANN and the per-entry placement keeps it aligned with the transcript it describes.
