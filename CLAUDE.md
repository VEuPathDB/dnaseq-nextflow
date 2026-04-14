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
nextflow run main.nf -entry loadSingleExperiment    -profile loadSingleExperiment
nextflow run main.nf -entry runTests                -profile tests
```

Docker is enabled by default in all profiles.

## Architecture

Three-tier structure: `main.nf` → `workflows/` → `modules/`

### Workflows

| Workflow | Purpose | Key modules |
|---|---|---|
| `processSingleExperiment` | Per-strain: FASTQ → consensus FASTA + VCF + coverage | preprocessing.nf, alignment.nf, snp.nf, cnv.nf |
| `mergeExperiments` | Multi-strain: merge VCFs, annotate variants, generate DB load files | mergeExperiments.nf |
| `loadSingleExperiment` | Load indel/ploidy/CNV data into GUS database | loadSingleExperiment.nf |
| `runTests` | Perl Test2::V0 test suite | runTests.nf |

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
  loadSingleExperiment.nf
modules/
  preprocessing.nf alignment.nf snp.nf cnv.nf
  mergeExperiments.nf loadSingleExperiment.nf runTests.nf
bin/
  processSequenceVariations.jl   # Core variation annotation (Julia)
  makeSnpFile.pl maskGenome.pl fixSeqId.pl
  calculatePloidy.pl calculateGeneCNVs.pl
  addFeatureIdsToVariation.pl addExtDbRlsIdToVariation.pl
testing/t/                       # Perl test files
testing/lib/                     # Test utilities
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
- `veupathdb/shortreadaligner:1.0.0` — BWA, samtools, Picard, GATK3, FreeBayes, bcftools, bedtools, Julia 1.10.10, Perl/BioPerl, SnpEff
- `veupathdb/dnaseqanalysis:1.0.0` — Trimmomatic, htseq-count

Julia deps (precompiled): `SQLite.jl`

## Testing

```bash
nextflow run main.nf -entry runTests -profile tests
```

Tests in `testing/t/` use Perl's `Test2::V0` framework, run via `prove`.
