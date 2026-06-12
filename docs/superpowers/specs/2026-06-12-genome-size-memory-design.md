# Genome-Size-Based Dynamic Memory — Design Spec

**Date:** 2026-06-12
**Branch:** parameterize_memory

## Problem

Memory limits for `bwaMem` and `runFreebayes` are currently calculated in the Perl workflow step (`MakeDnaSeqNextflowConfig.pm`) and written as hardcoded LSF `clusterOptions` MB values into the generated `nextflow.config`. `picard` and `gatk` have hardcoded `memory { 5.GB * task.attempt }` directives in their module files. All memory logic is scattered and the pipeline itself has no knowledge of its own resource requirements.

## Goal

Move all memory calculation logic into the Nextflow pipeline. Pass `genomeSize` (bytes) as a param from the Perl step. Use a single shared closure in `nextflow.config` to compute memory for all four affected processes.

## Constraints

- Input FASTQs are capped at 30× coverage, so BAM size scales predictably with genome size.
- LSF executor is used in production; Nextflow's `memory` directive maps to `-M` / `rusage[mem=]` automatically.
- `picard` and `gatk` scripts derive JVM heap from `task.memory` (`task.memory.toGiga() - 1`) — this must continue to work.

## Design

### 1. `MakeDnaSeqNextflowConfig.pm`

- Keep the genome size calculation (read FASTA line by line, sum sequence lengths).
- Emit `genomeSize = <bytes>` as a `params` entry in the generated config.
- Remove all memory math (power-of-2 rounding, `$bwaDefaultMemMb`, `$bwaRetryMemMb`).
- Remove `$bwaBlock` and `$freebayesBlock` string construction and their emission into the `process` block.
- Remove `maxMemoryGigs` param handling.

### 2. `nextflow.config`

Add `genomeSize` to the `params` block (default `0` — the Perl step always provides it).

Add a shared closure and four `withName` process blocks to the `process` block:

```groovy
def memForProcess = { genomeSizeBytes, multiplier, baseGb, attempt ->
    def genomeGb = genomeSizeBytes / 1_000_000_000
    def rawGb    = (genomeGb * multiplier) + baseGb
    def safeGb   = rawGb * 1.25
    def memGb    = 1; while (memGb < safeGb) memGb *= 2
    return "${memGb * attempt} GB"
}

process {
  executor = '...'
  queue    = '...'

  // To manually override memory for debugging, replace the closure with a fixed value, e.g.:
  //   withName: 'bwaMem' { memory = '32 GB' }
  withName: 'bwaMem'       { memory = { memForProcess(params.genomeSize, 3.3, 2, task.attempt) } }
  withName: 'runFreebayes' { memory = { memForProcess(params.genomeSize, 2.0, 2, task.attempt) } }
  withName: 'picard'       { memory = { memForProcess(params.genomeSize, 1.0, 2, task.attempt) } }
  withName: 'gatk'         { memory = { memForProcess(params.genomeSize, 1.0, 2, task.attempt) } }
}
```

### 3. `modules/alignment.nf`

- Remove `memory { 5.GB * task.attempt }` from `picard` and `gatk` process blocks.
- `jvmMem` derivation (`task.memory.toGiga() - 1`) stays unchanged — it reads from `task.memory` which is now set by the config closure.

### 4. Memory Formulas

All formulas apply: `round_up_to_power_of_2((genomeGb × multiplier + baseGb) × 1.25) × task.attempt`

| Process | Multiplier | Base GB | Rationale |
|---|---|---|---|
| `bwaMem` | 3.3 | 2 | SA index ≈ 3–4× genome in RAM |
| `runFreebayes` | 2.0 | 2 | Loads genome + 30× BAM regions |
| `picard` | 1.0 | 2 | Sort buffer ≈ BAM size ≈ genome at 30× coverage |
| `gatk` | 1.0 | 2 | Same memory profile as Picard |

Retries double the allocation via `task.attempt`.

### 5. Example Allocations

| Genome Size | bwaMem | FreeBayes | Picard | GATK |
|---|---|---|---|---|
| 5 MB (bacteria) | 4 GB | 4 GB | 4 GB | 4 GB |
| 100 MB (small euk) | 4 GB | 4 GB | 4 GB | 4 GB |
| 500 MB (mid euk) | 4 GB | 4 GB | 4 GB | 4 GB |
| 2 GB (large euk) | 16 GB | 8 GB | 4 GB | 4 GB |
| 3 GB (human-scale) | 16 GB | 8 GB | 8 GB | 8 GB |

## What Does Not Change

- LSF executor/queue configuration
- `errorStrategy` and `maxRetries` on `runFreebayes`
- `jvmMem` derivation in picard/gatk scripts
- All other pipeline logic

## Debugging Override

To pin a process to a fixed memory value in the generated config:
```groovy
withName: 'bwaMem' { memory = '32 GB' }
```
Replace the closure with a plain string — Nextflow accepts either form.
