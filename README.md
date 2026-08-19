# dnaseq-nextflow

A Nextflow DSL2 pipeline for DNA sequencing analysis: FASTQ reads to alignment, variant calling, CNV estimation, and multi-strain merge for downstream database loading.

## Overview

This pipeline processes raw whole-genome sequencing reads for individual strains/isolates and merges the results across strains into a form suitable for loading into VEuPathDB's GUS databases. It has two stages, exposed as separate entry points: `processSingleExperiment` aligns and calls variants for one strain at a time, producing a consensus FASTA, VCFs, indel and coverage tracks, ploidy, and gene-level CNV estimates; `mergeExperiments` combines the per-strain outputs of one or more `processSingleExperiment` runs, annotates variants against coding-sequence and indel databases, and runs SnpEff functional annotation to produce the variation/allele/product data files used by the VEuPathDB data-loading pipeline.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- [Docker](https://www.docker.com/) (enabled by default in every profile)

## Usage

```bash
# Per-strain analysis (default entry point)
nextflow run VEuPathDB/dnaseq-nextflow -r main -profile processSingleExperiment -resume -C <config>

# Per-strain analysis, named entry point
nextflow run VEuPathDB/dnaseq-nextflow -r main -entry processSingleExperiment -profile processSingleExperiment -resume -C <config>

# Multi-strain merge
nextflow run VEuPathDB/dnaseq-nextflow -r main -entry mergeExperiments -profile mergeExperiments -resume -C <config>
```

### Entry points

- **`processSingleExperiment`** (also the default workflow): Takes paired- or single-end FASTQs from an nf-core-style samplesheet and runs QC (FastQC, Trimmomatic) → alignment (BWA-MEM2, Picard dedup, GATK indel realignment) → variant calling (FreeBayes, consensus FASTA generation with coverage masking) → CNV/coverage analysis (genome coverage, htseq-count → TPM → ploidy and gene CNVs) → windowed SNP/heterozygous-SNP density and normalized coverage tracks.
- **`mergeExperiments`**: Combines per-strain indels, VCFs, and coverage BEDs from one or more `processSingleExperiment` runs; builds genomic-indel and coding-sequence/coding-indel SQLite databases from the consensus FASTAs and GTF; annotates variants with `bin/processSequenceVariations.jl`; and runs SnpEff for functional annotation. Outputs are intended for loading into a GUS/VEuPathDB database.

## Key parameters

### `processSingleExperiment`

| Parameter | Description |
|---|---|
| `samplesheet` | nf-core-format CSV (`sample`, `fastq_1`, `fastq_2`; `fastq_2` empty for single-end) |
| `genomeFastaFile` | Reference genome FASTA |
| `gtfFile` | Gene annotation GTF |
| `footprintFile` | Gene footprints file for CNV estimation |
| `geneSourceIdOrthologFile` | Gene source ID / ortholog mapping TSV |
| `chrsForCalcFile` | Chromosomes to include in ploidy calculation |
| `minCoverage` | Minimum depth required for variant calling and consensus masking |
| `ploidy` | Expected ploidy (heterozygous SNP tracks are skipped when `1`) |
| `winLen` | Window size (bp) for density and windowed coverage tracks |
| `bwaThreads` | Threads for BWA-MEM2 |
| `outputDir` | Output directory |

### `mergeExperiments`

| Parameter | Description |
|---|---|
| `relativeConsensusFilePattern` | Glob for per-strain `*_consensus.fa.gz` files |
| `vcfFiles` | Glob for per-strain VCF files |
| `indelsFiles` | Glob for per-strain `indels.tsv` files |
| `coverageFiles` | Glob for per-strain coverage BED files |
| `genomeFastaFile` | Reference genome FASTA |
| `gtfFile` | Gene annotation GTF |
| `cacheFile` / `vcfCacheFile` | VCF annotation cache from a previous run, to avoid re-annotating known variants |
| `undoneStrains` | Strains to exclude from annotation |
| `reference_strain` | Reference strain name |
| `outputDir` | Output directory |

### Global

| Parameter | Description |
|---|---|
| `benchmark` | When `true`, emits timing/call-count reports for the merge pipeline's Julia steps to `.command.err` |

## Output

**`processSingleExperiment`** produces, per strain: a coverage-masked consensus FASTA (`*_consensus.fa.gz`), the full FreeBayes VCF, an indel table (`indels.tsv`), per-position coverage BED, ploidy and gene-CNV estimates, coverage/normalized-coverage/SNP-density/het-SNP-density bigWig tracks, and a merged alignment-stats TSV.

**`mergeExperiments`** produces a merged, SnpEff-annotated multi-strain VCF along with variation, allele, and product DAT files formatted for loading into a GUS/VEuPathDB database.

## Containers

| Image | Used for |
|---|---|
| `veupathdb/shortreadaligner:1.0.1` | BWA-MEM2, samtools, Picard, GATK3, FreeBayes, bcftools, bedtools, Julia, Perl/BioPerl, SnpEff |
| `veupathdb/dnaseqanalysis:1.0.0` | Trimmomatic, htseq-count |

## Testing

Tests are Julia and Python unit tests under `testing/t/`, plus bash-driven characterization tests, run inside the `veupathdb/dnaseqanalysis` image (which carries Julia + SQLite.jl and Python + cyvcf2 + pytest):

```bash
docker run --rm --pull always -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:latest bash -c '
  for t in testing/t/*.jl; do julia "$t"; done   # Julia unit tests
  python3 -m pytest testing/t/                    # Python unit tests
  for t in testing/t/*.t.sh; do bash "$t"; done    # bash suites
'
```

The end-to-end tests in `testing/t/test_mergeExperiments_e2e.py` skip by default; run them against a completed `mergeExperiments` run with `python3 -m pytest testing/t/test_mergeExperiments_e2e.py --run-dir /path/to/run`. `testing/bin/compareMergeOutputs.py` compares two `mergeExperiments` output directories for equivalence modulo strain-id numbering.
