# Single-sample support in mergeExperiments

**Date**: 2026-08-01
**Status**: Approved, ready for implementation plan

## Problem

`mergeExperiments` assumes two or more samples. With a single experiment
contributing a single sample, two processes fail outright and a third code path
silently produces wrong input for everything downstream.

Verified inside `veupathdb/dnaseqanalysis:latest`:

| Command | Single-input behavior |
|---|---|
| `bcftools merge --merge both -O z -o merged.vcf.gz a.vcf.gz` | exit 1, usage dump (expects two or more files) |
| `bedtools unionbedg -names S1 -filler 0 -i c.bed` | exit 1, `Error: Only a single BedGraph file was specified. Nothing to combine, exiting.` |

`workflows/mergeExperiments.nf:29-31` already dodges the VCF crash:

```groovy
allVcfsBranch = allVcfs.branch { single: it.size() == 1; multiple: true }
mergedVcf     = allVcfsBranch.single.map { it[0] }
                  .mix(mergeVcfs(allVcfsBranch.multiple))
```

The `single` leg passes the **raw** per-sample VCF straight to `processSeqVars`,
skipping the two transformations `mergeVcfs` performs besides merging:

- `bcftools annotate -x "INFO,FORMAT/GL,FORMAT/DPR"` — strips freebayes INFO and
  the GL/DPR FORMAT fields
- `bcftools norm -m -any --multi-overlaps .` — splits multiallelic records into
  biallelic rows

So an n=1 run feeds `processSequenceVariations.jl` a structurally different VCF
than every n>1 run. This is live today, not a hypothetical.

`mergeCoverageBeds` has no workaround at all — an n=1 run dies there.

## Root cause

The output *contract* of each process (a normalized, INFO-stripped,
biallelic-split `merged.vcf.gz`; a `coverage.tsv` with a header plus per-sample
union intervals) is currently defined in two places: the process script, for the
n>1 case, and the workflow-level branch, for n=1. Two definitions of one contract
drift. They already have.

Input arity is an implementation detail of *how* a process satisfies its
contract. It does not belong in the workflow DAG.

## Design

### 1. `workflows/mergeExperiments.nf`

Delete the `branch`/`mix`. The workflow becomes arity-agnostic:

```groovy
mergedVcf   = mergeVcfs(vcfs_qch.collect())
coverageTsv = mergeCoverageBeds(coverages_qch.collect())
```

### 2. `modules/mergeExperiments.nf` — `mergeVcfs`

The normalization loop is unchanged: every input still gets `annotate -x` +
`norm -m -any --multi-overlaps .` + `index --tbi`. Only the final combine step
becomes arity-aware:

- **0 inputs** — write a clear error to stderr, exit non-zero.
- **1 input** — `mv` the single `*.norm.vcf.gz` to `merged.vcf.gz`.
- **≥2 inputs** — `bcftools merge --merge both` exactly as today.

Equivalence of the n=1 path was verified empirically: `bcftools merge` does
**not** re-collapse rows that `norm -m -any` split within a single file. A
multiallelic `A -> G,T` split into two biallelic rows stays two rows through the
merge. The n=1 output is therefore shape-identical to what merge would produce
for that strain, minus the extra `./.` columns contributed by other samples —
which by definition do not exist at n=1.

### 3. `modules/mergeExperiments.nf` — `mergeCoverageBeds`

Same structure. Sample-name derivation and header construction are unchanged
(`basename "$f" _coverage.bed.gz`), so the header is identical in every case.
Only the body differs:

- **0 inputs** — error to stderr, exit non-zero.
- **1 input** — `zcat` the single bed, appended after the header.
- **≥2 inputs** — `bedtools unionbedg` exactly as today.

The per-sample beds from `modules/snp.nf` are 4-column (`chrom start end mean`,
from `bedtools genomecov -bga | awk minCoverage | bedtools merge -c 4 -o mean`).
For a single file, `unionbedg` output is byte-identical to the file's own
contents, so `zcat` is exact — not an approximation.

## Strain ids are run-order dependent

Generated strain ids follow the order samples appear in the merged VCF header,
which follows the order the collected files were staged, which follows channel
order from `Channel.fromPath`. That order is not stable across runs. The existing
baseline at `~/dnaseq_test/merge/output.orig/` shows
`1=LV39, 2=Friedlin_resequence, 3=Seidman751, 4=lmajFriedlin` — neither
alphabetical nor otherwise reproducible.

Artifacts carrying ids:

| Artifact | Id-bearing content |
|---|---|
| `sample.dat` | `strain_id` → `sample_name` rows |
| `allele.dat` | `strain_ids` sets, e.g. `{1,2,4}` |
| `transcript_product.dat` | `downstream_of_frameshift_strain_ids` |
| `hsss_readFreq*/` | per-strain files *named* by index, plus `strainIdToName.dat` |
| `merged.ann.vcf.gz` | per-sample columns, in header order |

`variationFeature.dat` and `snpeff.dat` carry no ids.

Note the two numberings differ: `hsss_readFreq*/strainIdToName.dat` puts the
reference first (`1 lmajFriedlin`) while `sample.dat` puts it last (`4`). The
comparator must read each mapping from its own file rather than assume one
scheme.

**Deliberately deferred:** making this ordering deterministic (sorting samples by
name before merging) would make the pipeline's output reproducible and reduce
every future regression check to a plain `diff`. That was considered and
explicitly deferred — it is not part of this work. The consequence accepted here
is that comparing any two `mergeExperiments` runs requires the comparator below,
indefinitely.

## Testing

### Unit / characterization

New `testing/t/singleSampleMerge.t.sh`, run in `veupathdb/dnaseqanalysis:latest`
alongside the existing `*.t.sh` suites. Synthetic fixtures, no pipeline run.

Asserts:

1. n=1 `merged.vcf.gz` — INFO stripped, `FORMAT/GL` and `FORMAT/DPR` absent,
   multiallelic input emerges as separate biallelic rows, output is
   bgzip-compressed and tabix-indexable.
2. n=1 and n≥2 `coverage.tsv` share the same header shape
   (`chrom start end <name>...`) and column count `3 + nSamples`.
3. n=1 `coverage.tsv` data rows equal the input bed's rows.
4. n=0 — both script bodies exit non-zero with a message naming the missing
   input.

### Comparator

New `testing/bin/compareMergeOutputs.py`, taking two output directories and
reporting whether they are equivalent modulo strain-id permutation. Deliberately
placed outside `testing/t/` so pytest does not collect it — it is an instrument,
not a test. It reads each run's own id mappings and canonicalizes by sample name:

| Artifact | Comparison |
|---|---|
| `sample.dat` | compare as a set of sample names; ignore id assignment |
| `allele.dat` | rewrite `strain_ids` sets to name-sorted names, then compare |
| `transcript_product.dat` | same rewrite for `downstream_of_frameshift_strain_ids` |
| `variationFeature.dat`, `snpeff.dat` | compare directly (no ids) |
| `hsss_readFreq*/` | map filenames to strain names via each dir's `strainIdToName.dat`, compare contents pairwise by name; compare `contigIdToSourceId.dat` and `referenceGenome.dat` directly |
| `merged.ann.vcf.gz` | compare the record body only (`bcftools view -H`) with samples forced to a canonical order via `bcftools view -s` |

The VCF body-only comparison is required because bcftools and snpEff stamp
version and command-line lines into the `##` header, so the compressed files
differ even when the data does not.

Row order within a file is not guaranteed to be stable either (rows at one
position may be emitted in strain-iteration order), so the comparator sorts rows
before comparing.

### End-to-end

Two prepared run directories, invoked from within each:

```bash
nextflow -C $PWD/nextflow.config run \
  /home/jbrestel/workspaces/dataLoad/project_home/dnaseq-nextflow/main.nf \
  -entry mergeExperiments -profile mergeExperiments
```

- `~/dnaseq_test/single_sample` — one sample, `WhitePaper/Seidman751`.
- `~/dnaseq_test/merge` — four samples across `Mottran` and `WhitePaper`,
  including that same `Seidman751`. Baseline output already present at
  `output.orig/`.

**Multi-sample regression — must show zero change.** The n≥2 code path is
untouched by design. Rerun the `merge` directory and compare the new output
against `output.orig/` with the comparator. Any reported difference means the
n≥2 path was altered and the change is wrong.

**Single-sample checks:**

- The run completes through `parseSnpEffAnnotations`.
- `sample.dat` contains the one strain plus `reference_strain`
  (`lmajFriedlin`). `reference_strain` is external —
  `processSequenceVariations.jl` appends it to `all_strains` rather than
  selecting it from them — so n=1 is semantically well-formed: one strain plus
  the reference.
- `hsss_readFreq{20,40,60,80}` are produced.
- Seidman751's per-strain values are consistent between the two runs where they
  should be. Frequency and count columns legitimately differ (denominator is 1
  strain vs 4), so this is a sanity read, not a diff.

## Out of scope

- No changes to `bin/processSequenceVariations.jl`. If the e2e run shows the
  Julia code itself misbehaving at n=1 (degenerate major/minor allele calls,
  frequency edge cases), that gets reported and specced separately.
- Making strain-id assignment deterministic — see the deferred note above.
- No changes to `processSingleExperiment`.
- No changes to the other `mergeExperiments` processes (`makeGenomicIndelDb`,
  `makeCodingData`, `snpEff`, `parseSnpEffAnnotations`) — none of them require
  ≥2 inputs.
