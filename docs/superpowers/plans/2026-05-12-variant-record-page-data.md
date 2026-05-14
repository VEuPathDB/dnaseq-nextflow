# Variant Record Page Data Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the `mergeExperiments` pipeline to emit `sample.dat` and `snpeff.dat`, and add `strain_ids integer[]` to `allele.dat`, enabling per-strain allele breakdown and SnpEff impact display on the variant record page.

**Architecture:** `processSequenceVariations.jl` gains a `sample_id_map` (built from `all_strains` in the VCF header) and uses it to append a `strain_ids` array column to every `allele.dat` row; it also writes `sample.dat` once at startup. A new post-SnpEff Python script `parseSnpEffAnnotations.py` parses the `ANN` INFO field from `merged.ann.vcf.gz` and emits `snpeff.dat`. Both new artifacts are wired into the Nextflow `mergeExperiments` workflow.

**Tech Stack:** Julia 1.10, Python 3, Nextflow DSL2, bcftools (for VCF), PostgreSQL array syntax (`{1,2,3}`)

**Spec:** `docs/superpowers/specs/2026-05-12-variant-record-page-design.md`

---

## File Map

| Action | File | Responsibility |
|---|---|---|
| Modify | `bin/processSequenceVariations.jl` | Add `sample_id_map` to `ProcessingContext`, write `sample.dat`, add `strain_ids` to `allele.dat` |
| Create | `bin/parseSnpEffAnnotations.py` | Parse SnpEff ANN field → `snpeff.dat` |
| Modify | `modules/mergeExperiments.nf` | Add `sample.dat` output to `processSeqVars`; add `parseSnpEffAnnotations` process |
| Modify | `workflows/mergeExperiments.nf` | Wire `parseSnpEffAnnotations` after `snpEff` |
| Modify | `testing/t/handleVariantRecord.jl` | Tests for `write_allele_and_product_files` with `strain_ids` |
| Create | `testing/t/test_parseSnpEffAnnotations.py` | Tests for the new Python script |

---

## Task 1: Add `sample_id_map` to `ProcessingContext`

**Files:**
- Modify: `bin/processSequenceVariations.jl:408-418` (`ProcessingContext` struct)
- Modify: `bin/processSequenceVariations.jl:1035-1066` (`initialize_processing_context`)
- Test: `testing/t/handleVariantRecord.jl`

- [ ] **Step 1: Write the failing test**

Add to `testing/t/handleVariantRecord.jl` before the final `@testset` block:

```julia
@testset "initialize_processing_context builds sample_id_map from all_strains" begin
    # Use an empty/minimal args dict — we only care about sample_id_map here.
    # initialize_processing_context needs DB paths; use temp files for the SQLite dbs.
    tmp_dir = mktempdir()
    transcript_db_path = joinpath(tmp_dir, "transcripts.db")
    indel_db_path      = joinpath(tmp_dir, "indels.db")
    # Create minimal sqlite databases so the function doesn't error on open
    SQLite.DB(transcript_db_path) |> close
    SQLite.DB(indel_db_path)      |> close
    # Minimal GTF with no CDS entries
    gtf_path = joinpath(tmp_dir, "empty.gtf")
    write(gtf_path, "")

    args = Dict(
        "transcript_db"      => transcript_db_path,
        "indel_db"           => indel_db_path,
        "gtf_file"           => gtf_path,
        "reference_strain"   => "refStrain",
        "undone_strains_file"=> "",
    )
    all_strains = ["strainA", "strainB", "strainC"]
    ctx = initialize_processing_context(args, all_strains)

    @test ctx.sample_id_map["strainA"] == 1
    @test ctx.sample_id_map["strainB"] == 2
    @test ctx.sample_id_map["strainC"] == 3
    @test length(ctx.sample_id_map) == 3
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia testing/t/handleVariantRecord.jl 2>&1 | grep -E "FAIL|ERROR|sample_id_map"
```

Expected: ERROR — `type ProcessingContext has no field sample_id_map`

- [ ] **Step 3: Add `sample_id_map` field to `ProcessingContext`**

In `bin/processSequenceVariations.jl`, find the `ProcessingContext` struct (lines 408-418) and add the new field at the end:

```julia
# BEFORE (line 417):
    all_strains::Vector{String}
end

# AFTER:
    all_strains::Vector{String}
    sample_id_map::Dict{String,Int}
end
```

- [ ] **Step 4: Build `sample_id_map` in `initialize_processing_context`**

Find `initialize_processing_context` (around line 1037). The function ends with a constructor call. Add `sample_id_map` to the constructor. The function currently returns a `ProcessingContext(...)` call — add the field before the closing paren:

```julia
# Find the closing lines of initialize_processing_context, which look like:
        all_strains
    )
end

# Change to:
        all_strains,
        Dict{String,Int}(name => i for (i, name) in enumerate(all_strains))
    )
end
```

- [ ] **Step 5: Run test to verify it passes**

```bash
julia testing/t/handleVariantRecord.jl 2>&1 | grep -E "Test Summary|PASS|FAIL|ERROR"
```

Expected: all existing tests still pass, new testset shows `1 passed`.

- [ ] **Step 6: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "feat: add sample_id_map to ProcessingContext"
```

---

## Task 2: Write `sample.dat` from `main()`

**Files:**
- Modify: `bin/processSequenceVariations.jl` — add `write_sample_dat` helper and call it from `main()`
- Test: `testing/t/handleVariantRecord.jl`

- [ ] **Step 1: Write the failing test**

Add to `testing/t/handleVariantRecord.jl`:

```julia
@testset "write_sample_dat writes strain_id/sample_name TSV" begin
    tmp = tempname()
    write_sample_dat(["alpha", "beta", "gamma"], tmp)
    lines = readlines(tmp)
    @test lines[1] == "strain_id\tsample_name"
    @test lines[2] == "1\talpha"
    @test lines[3] == "2\tbeta"
    @test lines[4] == "3\tgamma"
    @test length(lines) == 4
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia testing/t/handleVariantRecord.jl 2>&1 | grep -E "FAIL|ERROR|write_sample_dat"
```

Expected: ERROR — `UndefVarError: write_sample_dat not defined`

- [ ] **Step 3: Add `write_sample_dat` function**

Add this function in `bin/processSequenceVariations.jl` after the `OutputWriters` struct definition (around line 426):

```julia
"""
    write_sample_dat(all_strains, path="sample.dat")

Write a two-column TSV (strain_id, sample_name) assigning sequential integer
IDs to samples in VCF column order.
"""
function write_sample_dat(all_strains::Vector{String}, path::String="sample.dat")
    open(path, "w") do fh
        write(fh, "strain_id\tsample_name\n")
        for (i, name) in enumerate(all_strains)
            write(fh, "$(i)\t$(name)\n")
        end
    end
end
```

- [ ] **Step 4: Call `write_sample_dat` from `main()`**

In `main()`, find line 1878:
```julia
    writers = open_output_writers(args["output_vcf"], args["output_cache"])
```

Add the call immediately after `all_strains` is obtained (line 1861). Insert after line 1862:
```julia
    write_sample_dat(all_strains)
```

- [ ] **Step 5: Run test to verify it passes**

```bash
julia testing/t/handleVariantRecord.jl 2>&1 | grep -E "Test Summary|PASS|FAIL|ERROR"
```

Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "feat: write sample.dat from all_strains VCF header"
```

---

## Task 3: Add `strain_ids` column to `allele.dat`

**Files:**
- Modify: `bin/processSequenceVariations.jl` — update `open_output_writers` header, `write_allele_and_product_files`, and its caller in `handle_variant_record!`
- Test: `testing/t/handleVariantRecord.jl`

- [ ] **Step 1: Write the failing test**

Add to `testing/t/handleVariantRecord.jl`:

```julia
@testset "write_allele_and_product_files emits strain_ids array in allele.dat" begin
    sample_id_map = Dict{String,Int}("s1" => 1, "s2" => 2, "s3" => 3)

    # Three hom variations: s1→A, s2→G, s3→G
    variations = [
        Variation("chr1", 100, "A", [], "s1", "A", "", "", "", 0, String[], ""),
        Variation("chr1", 100, "A", [], "s2", "G", "", "", "", 0, String[], ""),
        Variation("chr1", 100, "A", [], "s3", "G", "", "", "", 0, String[], ""),
    ]
    annotation = PositionAnnotation(0, "", 0, 0, 0, "", "")

    allele_buf  = IOBuffer()
    product_buf = IOBuffer()
    write_allele_and_product_files(allele_buf, product_buf, variations, annotation, "LmjF.01", 100, sample_id_map)

    lines = filter(!isempty, split(String(take!(allele_buf)), "\n"))
    @test length(lines) == 2   # one row per distinct allele

    # Find the A row and G row
    a_row = findfirst(l -> split(l, "\t")[3] == "A", lines)
    g_row = findfirst(l -> split(l, "\t")[3] == "G", lines)
    @test !isnothing(a_row)
    @test !isnothing(g_row)

    a_fields = split(lines[a_row], "\t")
    g_fields = split(lines[g_row], "\t")

    # strain_ids is the 8th column (index 8)
    @test a_fields[8] == "{1}"
    @test g_fields[8] == "{2,3}"   # sorted
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia testing/t/handleVariantRecord.jl 2>&1 | grep -E "FAIL|ERROR|strain_ids"
```

Expected: ERROR — wrong number of arguments to `write_allele_and_product_files`

- [ ] **Step 3: Update `allele.dat` header in `open_output_writers`**

Find line 1078:
```julia
    write(allele_fh,  "location\tseq_id\tallele\tdistinct_strain_count\tallele_count\tavg_coverage\tavg_percent\n")
```

Change to:
```julia
    write(allele_fh,  "location\tseq_id\tallele\tdistinct_strain_count\tallele_count\tavg_coverage\tavg_percent\tstrain_ids\n")
```

- [ ] **Step 4: Update `write_allele_and_product_files` signature and body**

Find the function at line 1333. Add `sample_id_map::Dict{String,Int}` parameter:

```julia
function write_allele_and_product_files(
    allele_fh::IO,
    product_fh::IO,
    variations::Vector{Variation},
    annotation::PositionAnnotation,
    seq_id::String,
    location::Int,
    sample_id_map::Dict{String,Int}
)
```

In the inner loop body (around line 1358), update the per-allele accumulation loop. Find:
```julia
    for (allele, entries) in allele_entries
        distinct_strains = Set{String}()
        sum_coverage = 0.0
        sum_percent  = 0.0

        for (strain, cov, pct) in entries
            push!(distinct_strains, strain)
            sum_coverage += cov
            sum_percent  += pct
        end

        n = length(entries)
        write(allele_fh, join([
            string(location),
            seq_id,
            allele,
            string(length(distinct_strains)),
            string(n),
            @sprintf("%.2f", sum_coverage / n),
            @sprintf("%.2f", sum_percent  / n)
        ], "\t"), "\n")
    end
```

Replace with:
```julia
    for (allele, entries) in allele_entries
        distinct_strains = Set{String}()
        strain_ids       = Set{Int}()
        sum_coverage = 0.0
        sum_percent  = 0.0

        for (strain, cov, pct) in entries
            push!(distinct_strains, strain)
            sid = get(sample_id_map, strain, 0)
            sid > 0 && push!(strain_ids, sid)
            sum_coverage += cov
            sum_percent  += pct
        end

        n        = length(entries)
        ids_str  = "{" * join(sort(collect(strain_ids)), ",") * "}"
        write(allele_fh, join([
            string(location),
            seq_id,
            allele,
            string(length(distinct_strains)),
            string(n),
            @sprintf("%.2f", sum_coverage / n),
            @sprintf("%.2f", sum_percent  / n),
            ids_str
        ], "\t"), "\n")
    end
```

- [ ] **Step 5: Update the caller in `handle_variant_record!`**

Find line 1804:
```julia
        write_allele_and_product_files(writers.allele_fh, writers.product_fh, all_vars, annotation, seq_id, location)
```

Change to:
```julia
        write_allele_and_product_files(writers.allele_fh, writers.product_fh, all_vars, annotation, seq_id, location, ctx.sample_id_map)
```

- [ ] **Step 6: Run tests to verify they pass**

```bash
julia testing/t/handleVariantRecord.jl 2>&1 | grep -E "Test Summary|PASS|FAIL|ERROR"
```

Expected: all tests pass including the new strain_ids testset.

- [ ] **Step 7: Commit**

```bash
git add bin/processSequenceVariations.jl testing/t/handleVariantRecord.jl
git commit -m "feat: add strain_ids array column to allele.dat"
```

---

## Task 4: Create `parseSnpEffAnnotations.py`

**Files:**
- Create: `bin/parseSnpEffAnnotations.py`
- Create: `testing/t/test_parseSnpEffAnnotations.py`

- [ ] **Step 1: Write the failing test**

Create `testing/t/test_parseSnpEffAnnotations.py`:

```python
"""Tests for bin/parseSnpEffAnnotations.py"""
import gzip
import os
import subprocess
import tempfile

import pytest

SCRIPT = os.path.join(os.path.dirname(__file__), "../../bin/parseSnpEffAnnotations.py")


def run_script(vcf_content: str, compressed: bool = True) -> list[str]:
    """Write VCF content to a temp file, run script, return output lines (no header)."""
    with tempfile.TemporaryDirectory() as tmp:
        vcf_path = os.path.join(tmp, "test.vcf.gz" if compressed else "test.vcf")
        out_path  = os.path.join(tmp, "snpeff.dat")
        if compressed:
            with gzip.open(vcf_path, 'wt') as f:
                f.write(vcf_content)
        else:
            with open(vcf_path, 'w') as f:
                f.write(vcf_content)
        subprocess.run(
            ["python3", SCRIPT, "--vcf", vcf_path, "--output", out_path],
            check=True
        )
        with open(out_path) as f:
            lines = f.read().strip().split('\n')
        return lines[1:]  # skip header


MINIMAL_HEADER = "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

ANN_MODERATE = "T|missense_variant|MODERATE|geneA|geneA_id|transcript|tx1|protein_coding|1/3|c.100A>T|p.Ser34Cys|100/900|100/900|34/300|.|"
ANN_HIGH     = "A|stop_gained|HIGH|geneA|geneA_id|transcript|tx1|protein_coding|1/3|c.50G>A|p.Trp17*|50/900|50/900|17/300|.|"


def test_parses_single_ann_entry():
    vcf = MINIMAL_HEADER + f"LmjF.01\t3745\t.\tC\tT\t.\t.\tANN={ANN_MODERATE}\n"
    rows = run_script(vcf)
    assert len(rows) == 1
    fields = rows[0].split('\t')
    assert fields[0] == "3745"           # location
    assert fields[1] == "LmjF.01"        # seq_id
    assert fields[2] == "T"              # allele
    assert fields[3] == "tx1"            # transcript_id
    assert fields[4] == "MODERATE"       # snpeff_impact
    assert fields[5] == "missense_variant"  # snpeff_effect


def test_parses_multiple_ann_entries_same_position():
    vcf = MINIMAL_HEADER + f"LmjF.01\t3745\t.\tG\tT,A\t.\t.\tANN={ANN_MODERATE},{ANN_HIGH}\n"
    rows = run_script(vcf)
    assert len(rows) == 2
    alleles = {r.split('\t')[2] for r in rows}
    impacts = {r.split('\t')[4] for r in rows}
    assert alleles == {"T", "A"}
    assert impacts == {"MODERATE", "HIGH"}


def test_skips_entries_with_empty_transcript():
    no_tx = "T|missense_variant|MODERATE|geneA|geneA_id|transcript||protein_coding|1/3|.|.|.|.|.|.|"
    vcf   = MINIMAL_HEADER + f"LmjF.01\t100\t.\tC\tT\t.\t.\tANN={no_tx}\n"
    rows  = run_script(vcf)
    assert rows == [""]   # no data rows (empty string from strip+split on empty file body)


def test_deduplicates_same_allele_transcript_pair():
    dup = f"{ANN_MODERATE},{ANN_MODERATE}"
    vcf = MINIMAL_HEADER + f"LmjF.01\t200\t.\tC\tT\t.\t.\tANN={dup}\n"
    rows = [r for r in run_script(vcf) if r]
    assert len(rows) == 1


def test_handles_uncompressed_vcf():
    vcf  = MINIMAL_HEADER + f"LmjF.01\t500\t.\tC\tT\t.\t.\tANN={ANN_MODERATE}\n"
    rows = run_script(vcf, compressed=False)
    assert len(rows) == 1
    assert rows[0].split('\t')[4] == "MODERATE"


def test_skips_records_without_ann():
    vcf = MINIMAL_HEADER + "LmjF.01\t999\t.\tC\tT\t.\t.\tDP=30\n"
    rows = [r for r in run_script(vcf) if r]
    assert rows == []
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd testing/t && python3 -m pytest test_parseSnpEffAnnotations.py -v 2>&1 | head -20
```

Expected: all FAIL — `FileNotFoundError` or `subprocess.CalledProcessError` (script doesn't exist yet).

- [ ] **Step 3: Write `bin/parseSnpEffAnnotations.py`**

```python
#!/usr/bin/env python3
"""Parse SnpEff ANN INFO field from annotated VCF and write snpeff.dat.

Output columns: location, seq_id, allele, transcript_id, snpeff_impact, snpeff_effect
One row per unique (location, seq_id, allele, transcript_id) combination.
"""
import argparse
import gzip


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf",    required=True, help="SnpEff-annotated VCF (.vcf or .vcf.gz)")
    p.add_argument("--output", required=True, help="Output path for snpeff.dat")
    return p.parse_args()


def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    args = parse_args()
    seen = set()

    with open_vcf(args.vcf) as vcf_fh, open(args.output, "w") as out:
        out.write("location\tseq_id\tallele\ttranscript_id\tsnpeff_impact\tsnpeff_effect\n")
        for line in vcf_fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            seq_id   = fields[0]
            location = fields[1]
            info     = fields[7]

            ann_value = None
            for token in info.split(";"):
                if token.startswith("ANN="):
                    ann_value = token[4:]
                    break
            if ann_value is None:
                continue

            for entry in ann_value.split(","):
                parts = entry.split("|")
                if len(parts) < 7:
                    continue
                allele        = parts[0]
                effect        = parts[1]
                impact        = parts[2]
                transcript_id = parts[6]
                if not transcript_id or not impact or not effect:
                    continue
                key = (location, seq_id, allele, transcript_id)
                if key in seen:
                    continue
                seen.add(key)
                out.write(f"{location}\t{seq_id}\t{allele}\t{transcript_id}\t{impact}\t{effect}\n")


if __name__ == "__main__":
    main()
```

Make it executable:
```bash
chmod +x bin/parseSnpEffAnnotations.py
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd testing/t && python3 -m pytest test_parseSnpEffAnnotations.py -v 2>&1 | grep -E "PASSED|FAILED|ERROR"
```

Expected: all 6 tests PASSED.

- [ ] **Step 5: Commit**

```bash
git add bin/parseSnpEffAnnotations.py testing/t/test_parseSnpEffAnnotations.py
git commit -m "feat: add parseSnpEffAnnotations.py to emit snpeff.dat"
```

---

## Task 5: Nextflow wiring

**Files:**
- Modify: `modules/mergeExperiments.nf` — add `sample.dat` output to `processSeqVars`; add `parseSnpEffAnnotations` process
- Modify: `workflows/mergeExperiments.nf` — wire `parseSnpEffAnnotations` after `snpEff`

- [ ] **Step 1: Add `sample.dat` output to the `processSeqVars` process**

In `modules/mergeExperiments.nf`, find the `processSeqVars` output block (around line 162):
```groovy
  output:
    path 'output.vcf.gz', emit: outputVcf
    path 'output.vcf.gz.tbi', emit: outputVcfIndex
    path 'output_cache.tsv', emit: outputCache
    path 'variationFeature.dat', emit: variationFile
    path 'allele.dat', emit: alleleFile
    path 'product.dat', emit: productFile
```

Add `sample.dat`:
```groovy
  output:
    path 'output.vcf.gz', emit: outputVcf
    path 'output.vcf.gz.tbi', emit: outputVcfIndex
    path 'output_cache.tsv', emit: outputCache
    path 'variationFeature.dat', emit: variationFile
    path 'allele.dat', emit: alleleFile
    path 'product.dat', emit: productFile
    path 'sample.dat', emit: sampleFile
```

Also add a `publishDir` for `sample.dat` alongside the others (around line 147):
```groovy
  publishDir "$params.outputDir", mode: "copy", pattern: 'sample.dat'
```

Also add to the stub block:
```groovy
    touch sample.dat
```

- [ ] **Step 2: Add `parseSnpEffAnnotations` process to `modules/mergeExperiments.nf`**

Append after the `snpEff` process (after line 234):

```groovy
process parseSnpEffAnnotations {
  container 'veupathdb/shortreadaligner:1.0.0'
  publishDir "$params.outputDir", mode: "copy", pattern: 'snpeff.dat'

  input:
    path annVcf

  output:
    path 'snpeff.dat', emit: snpeffFile

  script:
    """
    parseSnpEffAnnotations.py --vcf $annVcf --output snpeff.dat
    """

  stub:
    """
    touch snpeff.dat
    """
}
```

- [ ] **Step 3: Add the include and wire in `workflows/mergeExperiments.nf`**

At the top of `workflows/mergeExperiments.nf`, add the import:
```groovy
include { parseSnpEffAnnotations } from '../modules/mergeExperiments.nf'
```

In the `main:` block of the `me` workflow, find the `snpEff` call:
```groovy
    snpEff(processSeqVarsResults.outputVcf, params.gtfFile, params.genomeFastaFile)
```

Add the call immediately after:
```groovy
    parseSnpEffAnnotations(snpEff.out[0])
```

(`snpEff` emits `merged.ann.vcf.gz` as its first unnamed output.)

- [ ] **Step 4: Verify the workflow stub-runs without error**

```bash
nextflow run main.nf -entry mergeExperiments -profile mergeExperiments -stub 2>&1 | tail -20
```

Expected: workflow completes with `N processes cached/completed`, no ERROR lines.

- [ ] **Step 5: Verify `sample.dat` and `snpeff.dat` appear in output dir after a stub run**

```bash
ls output/sample.dat output/snpeff.dat
```

Expected: both files present (empty from stub).

- [ ] **Step 6: Commit**

```bash
git add modules/mergeExperiments.nf workflows/mergeExperiments.nf
git commit -m "feat: wire sample.dat and snpeff.dat into mergeExperiments workflow"
```

---

## Self-Review Checklist

- [x] `sample.dat` — spec says "generated from samplesheet at load time"; implemented by deriving from VCF header `all_strains` (equivalent: the VCF already contains exactly the samples processed). No spec gap.
- [x] `allele.dat strain_ids` — covered by Task 3; PostgreSQL array syntax `{1,2,3}`.
- [x] `snpeff.dat` with impact + effect per transcript — covered by Task 4.
- [x] Nextflow wiring for all new files — covered by Task 5.
- [x] No TBDs or placeholders in any task.
- [x] Type consistency: `Dict{String,Int}` used for `sample_id_map` throughout; `strain_ids` written as string `"{1,2,3}"` (PostgreSQL COPY compatible).
- [x] Scope: page layout itself is a webapp concern, not implemented here — correctly out of scope.
