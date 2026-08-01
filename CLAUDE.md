# CLAUDE.md

Nextflow DSL2 pipeline for DNA sequencing analysis: FASTQ → alignment → variant calling → CNV → multi-strain merge → GUS database loading.

**Status**: Under construction, not used in production.

## Running Workflows

```bash
# Default (processSingleExperiment)
nextflow run main.nf -profile processSingleExperiment

# Named entry points
nextflow run main.nf -entry processSingleExperiment -profile processSingleExperiment
nextflow run main.nf -entry mergeExperiments        -profile mergeExperiments
```

Tests are not a Nextflow workflow — see [Testing](#testing).

Docker is enabled by default in all profiles.

## Architecture

Three-tier structure: `main.nf` → `workflows/` → `modules/`

### Workflows

| Workflow | Purpose | Key modules |
|---|---|---|
| `processSingleExperiment` | Per-strain: FASTQ → consensus FASTA + VCF + coverage | preprocessing.nf, alignment.nf, snp.nf, cnv.nf |
| `mergeExperiments` | Multi-strain: merge VCFs, annotate variants, generate DB load files | mergeExperiments.nf |

### processSingleExperiment stages
1. QC: FastQC, Trimmomatic
2. Alignment: BWA-MEM → Picard dedup → GATK indel realignment
3. Variant calling: FreeBayes → indel TSV → consensus + masked genome
4. CNV: bedtools coverage → htseq-count → TPM → ploidy + gene CNV
5. Windowed: SNP density, heterozygous SNP density, normalized coverage BigWigs

### mergeExperiments stages
1. Merge VCFs across strains (bcftools)
2. `bin/processSequenceVariations.jl` — annotates variants via SQLite transcript/indel DBs; outputs cache + variation/allele/product DAT files
3. Add GUS feature IDs, generate DB load files
4. SnpEff functional annotation

## Key Files

```
main.nf                          # Entry point, samplesheet parsing
nextflow.config                  # All profiles and parameters
workflows/
  processSingleExperiment.nf
  mergeExperiments.nf
modules/
  preprocessing.nf alignment.nf snp.nf cnv.nf
  mergeExperiments.nf
bin/
  processSequenceVariations.jl   # Core variation annotation (Julia)
  findValues.pl                  # Indel TSV extraction (snp.nf)
  calculatePloidy.pl calculateGeneCNVs.pl makeTpmFromHtseqCountsCNV.pl  # CNV (cnv.nf)
testing/t/                       # Julia (.jl) and Python (.py) tests
```

## Configuration

Key parameters in `nextflow.config` (profile-scoped):

| Parameter | Description |
|---|---|
| `samplesheet` | nf-core CSV (sample, fastq_1, fastq_2) |
| `genomeFastaFile` | Reference genome FASTA |
| `gtfFile` | Gene annotation GTF |
| `footprintFile` | Gene footprints for CNV |
| `minCoverage` | Min coverage for variant calling/masking |
| `ploidy` | Expected ploidy |

## Containers

Each process declares its own Docker image. Key images:
- `veupathdb/shortreadaligner:1.0.1` — BWA, samtools, Picard, GATK3, FreeBayes, bcftools, bedtools, Julia 1.10.10, Perl/BioPerl, SnpEff
- `veupathdb/dnaseqanalysis:1.0.0` — Trimmomatic, htseq-count

Julia deps (precompiled): `SQLite.jl`

## Testing

Tests are **not** run through Nextflow. They are Julia and Python unit tests in
`testing/t/`, plus a bcftools characterization test, run inside the
`veupathdb/dnaseqanalysis` image (which carries Julia + SQLite.jl, Python +
cyvcf2 + pytest, and bcftools). Use the `:latest` tag — Jenkins rebuilds it from
this repo's `Dockerfile` on every commit to `main`, so it stays in sync with the
code. `--pull always` refreshes any stale local copy:

```bash
docker run --rm --pull always -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest bash -c '
  for t in testing/t/*.jl; do julia "$t"; done   # Julia unit tests
  python3 -m pytest testing/t/                    # Python unit tests
  for t in testing/t/*.t.sh; do bash "$t"; done    # bash suites (bcftools merge contract, read-name normalization)
'
```

Individual suites can be run the same way (e.g. `julia testing/t/handleVariantRecord.jl`
or `python3 -m pytest testing/t/test_parseSnpEffAnnotations.py -v`).

**If your branch edits the `Dockerfile`**, `:latest` won't reflect it until the
branch merges to `main` and Jenkins rebuilds. Build locally and run against that
instead: `docker build -t dnaseqanalysis:dev .` then use `dnaseqanalysis:dev` above.

The 96 tests in `test_mergeExperiments_e2e.py` are end-to-end and **skip by
default** — they assert against the output files of a real `mergeExperiments`
run. To exercise them, point at a completed run:

```bash
python3 -m pytest testing/t/test_mergeExperiments_e2e.py --run-dir /path/to/nextflow/run
```

`testing/bin/compareMergeOutputs.py DIR_A DIR_B` compares two `mergeExperiments`
output directories for equivalence modulo strain-id permutation. Strain ids are
assigned in channel-staging order, so two runs over identical inputs can produce
identical data under different numbering — a plain `diff` reports differences
that are pure renumbering. Use this instead when checking that a change left
multi-sample output untouched. It canonicalizes the id sets in `allele.dat` and
`transcript_product.dat` by sample name, maps `hsss_readFreq*` filenames through
each directory's own `strainIdToName.dat` (a *different* numbering from
`sample.dat` — the reference strain comes first there), and compares the
annotated VCF's record body with sample columns reordered by name.
