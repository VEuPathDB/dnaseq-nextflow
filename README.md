# dnaseq-nextflow

> **Under construction — not used in production.**

Nextflow DSL2 pipeline for DNA sequencing analysis. Processes per-strain FASTQ files through alignment, variant calling, CNV, and coverage analysis (`processSingleExperiment`), then merges results across strains for downstream database loading (`mergeExperiments`).

---

## Quick Start

```bash
# Per-strain analysis (default)
nextflow run main.nf -profile processSingleExperiment

# Multi-strain merge
nextflow run main.nf -entry mergeExperiments -profile mergeExperiments
```

Tests are not run through Nextflow — see [Testing](#testing).

Docker is enabled by default in all profiles.

---

## processSingleExperiment

Runs per-strain: takes paired or single-end FASTQs from an nf-core samplesheet and produces a consensus FASTA, VCFs, indel TSV, coverage tracks, ploidy, and gene CNV estimates.

```
┌──────────────────────────────────────────────────────────┐
│                    QC & Preprocessing                    │
│                                                          │
│   FASTQ ──► FastQC ──► FastQC check ──► Trimmomatic      │
└────────────────────────────┬─────────────────────────────┘
                             │
┌────────────────────────────▼─────────────────────────────┐
│                        Alignment                         │
│                                                          │
│   BWA-MEM ──► Picard (dedup) ──► GATK (indel realign)    │
└──────────────────┬──────────────────────┬────────────────┘
                   │                      │
┌──────────────────▼──────┐  ┌────────────▼───────────────┐
│     Variant Calling     │  │       CNV / Coverage        │
│                         │  │                             │
│  FreeBayes              │  │  genomecov                  │
│  filterAndSplitVcf      │  │    └──► coverage bigwig     │
│    ├── SNPs VCF         │  │                             │
│    ├── Indels VCF       │  │  htseqCount ──► TPM         │
│    └── Consensus VCF    │  │    └──► ploidy & gene CNVs  │
│  makeIndelTSV           │  │                             │
│  makeCoverageBed        │  │  bedtoolsWindowed            │
│  consensus FASTA        │  │    └──► norm. cov. bigwig   │
│    (coverage-masked)    │  │                             │
└──────────────┬──────────┘  └─────────────────────────────┘
               │
┌──────────────▼──────────────────────────────────────────┐
│                    Density Tracks                        │
│                                                          │
│  SNP density bigwigs                                     │
│  Het SNP density bigwigs  (ploidy > 1 only)              │
└──────────────────────────────────────────────────────────┘

  Alignment stats (samtools + bedtools genomecov) ──► merged TSV
```

### Outputs

| File | Description |
|---|---|
| `*_consensus.fa.gz` | Per-strain consensus FASTA, low-coverage positions masked |
| `result.vcf.gz` | Full FreeBayes VCF (complex variants) |
| `indels.tsv` | Indel table (homozygous only; het indels excluded) |
| `coverage.bed.gz` | Per-position coverage BED |
| `*_Ploidy.txt` | Estimated ploidy per strain |
| `*_geneCNVs.txt` | Gene-level copy number estimates |
| `*.bw` | Coverage, normalised coverage, SNP density, het SNP density bigwigs |
| `alignment_stats.tsv` | Merged samtools + bedtools coverage stats across all samples |

### Parameters

| Parameter | Description |
|---|---|
| `samplesheet` | nf-core CSV (`sample`, `fastq_1`, `fastq_2`) |
| `genomeFastaFile` | Reference genome FASTA |
| `gtfFile` | Gene annotation GTF |
| `footprintFile` | Gene footprints file for CNV |
| `geneSourceIdOrthologFile` | Gene source ID / ortholog mapping TSV |
| `chrsForCalcFile` | Chromosomes to include in ploidy calculation |
| `minCoverage` | Minimum depth for variant calling and consensus masking |
| `ploidy` | Expected ploidy (het SNP tracks skipped when `1`) |
| `winLen` | Window size (bp) for density and windowed coverage tracks |
| `bwaThreads` | Threads for BWA-MEM |
| `outputDir` | Output directory |

---

## mergeExperiments

Takes the per-strain outputs from one or more `processSingleExperiment` runs and merges them across strains. Annotates variants against SQLite coding-sequence and indel databases, then runs SnpEff for functional annotation. Outputs are intended for loading into a GUS/VEuPathDB database.

### Steps

1. **Combine indels** — collect per-strain `indels.tsv` files and build a genomic indel SQLite database (`makeGenomicIndelDb`)
2. **Merge VCFs** — bcftools merge across all per-strain VCFs; single-strain inputs skip the merge step
3. **Merge coverage** — concatenate per-strain coverage BED files into a single TSV (`mergeCoverageBeds`)
4. **Build coding data** — derive coding-sequence and coding-indel SQLite databases from consensus FASTAs + GTF + reference genome (`makeCodingData`)
5. **Annotate variants** — `processSeqVars` runs `bin/processSequenceVariations.jl` against the SQLite DBs; produces annotated VCF and variation/allele/product DAT files
6. **SnpEff** — functional annotation of the merged VCF

### Inputs (from processSingleExperiment)

| Parameter | Description |
|---|---|
| `relativeConsensusFilePattern` | Glob for per-strain `*_consensus.fa.gz` files |
| `vcfFiles` | Glob for per-strain `result.vcf.gz` files |
| `indelsFiles` | Glob for per-strain `indels.tsv` files |
| `coverageFiles` | Glob for per-strain `*.coverage.bed.gz` files |

### Additional Parameters

| Parameter | Description |
|---|---|
| `genomeFastaFile` | Reference genome FASTA |
| `gtfFile` | Gene annotation GTF |
| `vcfCacheFile` | VCF cache from a previous run (avoid re-annotating known variants) |
| `undoneStrains` | Strains to exclude from annotation |
| `reference_strain` | Reference strain name |
| `outputDir` | Output directory |

---

## Containers

| Image | Used by |
|---|---|
| `veupathdb/shortreadaligner:1.0.0` | BWA, samtools, Picard, GATK3, FreeBayes, bcftools, bedtools, Julia 1.10, Perl/BioPerl, SnpEff |
| `veupathdb/dnaseqanalysis:1.0.0` | Trimmomatic, htseq-count |

---

## Testing

Tests are Julia and Python unit tests in `testing/t/`, plus a bcftools
characterization test. Run them inside the `veupathdb/dnaseqanalysis` image,
which carries Julia + SQLite.jl, Python + cyvcf2 + pytest, and bcftools:

```bash
docker run --rm -v "$PWD":/work -w /work veupathdb/dnaseqanalysis:1.1.1 bash -c '
  for t in testing/t/*.jl; do julia "$t"; done   # Julia unit tests
  python3 -m pytest testing/t/                    # Python unit tests
  bash testing/t/mergeBoth.t.sh                   # bcftools merge -m both contract
'
```

`pytest` is provided by the image as of the `1.1.1` rebuild; older pulls need
`pip3 install --break-system-packages pytest` first.
