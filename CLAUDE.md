# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a Nextflow-based DNA sequencing analysis pipeline for processing genomic variation data. The pipeline processes raw sequencing reads (FASTQ) through alignment, variant calling, and CNV analysis, then merges results across multiple strains for downstream database loading and analysis.

**Status**: Under construction, not used in production.

## Running Workflows

### Run the default workflow (processSingleExperiment)
```bash
nextflow run main.nf -profile processSingleExperiment
```

### Run specific workflows
```bash
# Process single experiment (alignment, variant calling, CNV)
nextflow run main.nf -entry processSingleExperiment -profile processSingleExperiment

# Merge results across multiple experiments
nextflow run main.nf -entry mergeExperiments -profile mergeExperiments

# Load single experiment results to database
nextflow run main.nf -entry loadSingleExperiment -profile loadSingleExperiment

# Run test suite
nextflow run main.nf -entry runTests -profile tests
```

### Docker execution
All workflows are designed to run in Docker containers. The profile configurations enable Docker by default (see `nextflow.config`).

## Architecture

### Workflow Organization

The codebase follows a three-tier structure:

1. **main.nf** - Entry point defining four named workflows: `processSingleExperiment`, `mergeExperiments`, `loadSingleExperiment`, and `runTests`
2. **workflows/** - High-level workflow orchestration that composes modules
3. **modules/** - Reusable process definitions grouped by function

### Three Primary Workflows

#### 1. processSingleExperiment (ps)
**Purpose**: Per-strain analysis from raw reads to variant calls and CNV data
**Input**: FASTQ files for individual strains (via nf-core samplesheet format)
**Output**: Consensus FASTA, VCF files, indel tables, ploidy estimates, gene CNVs, coverage bigWigs

**Pipeline stages**:
- Preprocessing: Quality control (FastQC), trimming (Trimmomatic)
- Alignment: BWA-MEM alignment, Picard deduplication, GATK realignment
- Variant calling: FreeBayes → SNP/indel separation → consensus genome generation
- CNV analysis: Coverage calculation, gene copy number estimation, ploidy determination
- Windowed analysis: SNP density, heterozygous SNP density, normalized coverage

**Key modules**: `preprocessing.nf`, `alignment.nf`, `snp.nf`, `cnv.nf`

#### 2. mergeExperiments (me)
**Purpose**: Combine multi-strain outputs and prepare for database loading
**Input**: Consensus FASTAs and VCF files from processSingleExperiment
**Output**: Merged VCF, annotated variation files, database load files, SnpEff annotations

**Pipeline stages**:
- Merge VCFs across all strains
- Process sequence variations using `bin/processSequenceVariations.jl` (Julia implementation replacing legacy Perl)
- Annotate variants with transcript/gene features
- Generate database load files (variation, product, allele tables)
- Run SnpEff for functional annotation

**Key modules**: `mergeExperiments.nf`

#### 3. loadSingleExperiment (ls)
**Purpose**: Load per-strain indel and CNV data into GUS database
**Input**: Indel TSV files, ploidy files, gene CNV files
**Key modules**: `loadSingleExperiment.nf`

### Module Structure

Modules are organized by analysis stage:
- **preprocessing.nf**: QC and trimming
- **alignment.nf**: Read alignment and BAM processing
- **snp.nf**: Variant calling and consensus generation
- **cnv.nf**: Copy number variation and coverage analysis
- **mergeExperiments.nf**: Multi-strain merging and annotation
- **loadSingleExperiment.nf**: Database loading
- **runTests.nf**: Test execution

### Key Processing Scripts

The `bin/` directory contains Perl and Julia scripts used by processes:

- **processSequenceVariations.jl**: Core variation annotation script (Julia rewrite, replaces processSequenceVariationsNew.pl)
  - Merges SNP file with cache file
  - Annotates coding variants with codon/product information via SQLite
  - Uses transcript and indel databases
  - Outputs: cache, snpFeature.dat, allele.dat, product.dat

- **Variant processing**: maskGenome.pl, makeSnpFile.pl, fixSeqId.pl
- **CNV calculation**: calculatePloidy.pl, calculateGeneCNVs.pl
- **Database utilities**: addFeatureIdsToVariation.pl, addExtDbRlsIdToVariation.pl

### Data Flow

```
FASTQ files (via samplesheet)
    ↓ (processSingleExperiment)
Per-strain: consensus FASTA + VCF + coverage
    ↓ (mergeExperiments)
Merged VCF + annotated variations + database files
    ↓ (loadSingleExperiment or database loading)
Populated GUS database
```

### Configuration

All parameters are defined in `nextflow.config` under profile-specific sections:
- Input/output directories
- Tool parameters (coverage thresholds, ploidy, variant calling parameters)
- Reference files (genome FASTA, GTF, footprints)
- Database connection details (for merge/load workflows)

Key parameters:
- `samplesheet`: Path to nf-core format CSV samplesheet (sample, fastq_1, fastq_2 columns)
- `minCoverage`: Minimum coverage threshold for variant calling and masking
- `ploidy`: Expected ploidy level
- `freebayesMinAltFraction`: Minimum allele frequency for variant calls

## Development

### Container and Dependencies

The Docker image (`veupathdb/shortreadaligner:1.0.0`) includes:
- Alignment tools: BWA, samtools, Picard, GATK
- Variant callers: FreeBayes, bcftools
- Analysis tools: bedtools, bedGraphToBigWig, htseq-count
- Languages: Perl (with BioPerl), Julia 1.10.10, Python
- VEuPathDB GUS framework components (for database loading)
- SnpEff for variant annotation

Julia dependencies (precompiled in image): SQLite.jl

### Testing

Tests are located in `testing/t/` and use Perl's Test2::V0 framework:
```bash
nextflow run main.nf -entry runTests -profile tests
```

Test utilities are in `testing/lib/`.

### Recent Refactoring

The Julia implementation (`bin/processSequenceVariations.jl`) was recently refactored to break up a 512-line main() function into modular functions. The variant calling has also been migrated from Varscan to FreeBayes.

## Input Data Requirements

### processSingleExperiment

**Samplesheet** (CSV format, nf-core standard):
- `sample`: Sample identifier (required, no spaces)
- `fastq_1`: Path to R1/forward reads file (required)
- `fastq_2`: Path to R2/reverse reads file (optional - leave empty for single-end)

Example samplesheet.csv:
```csv
sample,fastq_1,fastq_2
7G8,/path/to/7G8_R1.fastq.gz,/path/to/7G8_R2.fastq.gz
CS2,/path/to/CS2_R1.fastq.gz,/path/to/CS2_R2.fastq.gz
5.1,/path/to/5.1_SE.fastq.gz,
```

**Other required files**:
- Reference genome FASTA
- Gene annotation GTF file
- Gene footprints file
- Trimmomatic adapters file (optional, defaults to built-in adapters)

### mergeExperiments
- Consensus FASTA files (*.fa.gz) from processSingleExperiment
- VCF files (result.vcf.gz) from processSingleExperiment
- Coverage files (*.coverage.txt)
- Transcript SQLite database
- Indel SQLite database
- Cache file, undoneStrains file, gusConfig file
