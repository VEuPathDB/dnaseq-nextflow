# Genome-Size-Based Dynamic Memory Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move all memory limit logic for `bwaMem`, `runFreebayes`, `picard`, and `gatk` into the Nextflow pipeline using a genome-size-based closure, removing it from `MakeDnaSeqNextflowConfig.pm`.

**Architecture:** A shared `memForProcess` closure in `nextflow.config` computes memory as a function of `params.genomeSize` (bytes) and `task.attempt`. The Perl workflow step calculates genome size from the FASTA and passes it as a param; all memory math is removed from Perl. The `picard` and `gatk` module files drop their hardcoded `memory` directives since the config now provides them.

**Tech Stack:** Nextflow DSL2 (Groovy config closures), Perl, LSF cluster executor

---

## File Map

| File | Change |
|---|---|
| `nextflow.config` | Add global `params.genomeSize` default, `memForProcess` closure, and four `withName` process blocks |
| `modules/alignment.nf` | Remove `memory { 5.GB * task.attempt }` from `picard` (line 107) and `gatk` (line 147) |
| `ApiCommonWorkflow/Main/lib/perl/WorkflowSteps/MakeDnaSeqNextflowConfig.pm` | Keep genome size calculation, emit `genomeSize` param, remove `$bwaBlock`/`$freebayesBlock` memory logic and `maxMemoryGigs` handling |

---

### Task 1: Add `memForProcess` closure and process blocks to `nextflow.config`

**Files:**
- Modify: `nextflow.config` (top of file, before `profiles {`)

The current `nextflow.config` has only a `profiles { ... }` block. We add a global `params` default, a Groovy closure, and a top-level `process` block. Nextflow merges top-level and profile-scoped config at runtime.

- [ ] **Step 1: Open `nextflow.config` and prepend the following block before `profiles {`**

```groovy
// Default genomeSize — overridden by the cluster-side generated config.
// With genomeSize = 0 all processes floor to 4 GB, which is fine for local dev.
params {
  genomeSize = 0
}

// Compute memory as a function of genome size and retry attempt.
// Formula: round_up_to_power_of_2((genomeGb * multiplier + baseGb) * 1.25) * attempt
// With 30x coverage cap, BAM size scales with genome size, making this formula
// valid for all four processes below.
def memForProcess = { genomeSizeBytes, multiplier, baseGb, attempt ->
    def genomeGb = genomeSizeBytes / 1_000_000_000
    def rawGb    = (genomeGb * multiplier) + baseGb
    def safeGb   = rawGb * 1.25
    def memGb    = 1
    while (memGb < safeGb) memGb *= 2
    return "${memGb * attempt} GB"
}

process {
  // To manually override memory for debugging, replace the closure with a fixed value, e.g.:
  //   withName: 'bwaMem' { memory = '32 GB' }
  withName: 'bwaMem' {
    memory        = { memForProcess(params.genomeSize, 3.3, 2, task.attempt) }
    maxRetries    = 1
    errorStrategy = { task.exitStatus in 130..140 ? 'retry' : 'finish' }
  }
  withName: 'runFreebayes' {
    memory        = { memForProcess(params.genomeSize, 2.0, 2, task.attempt) }
    maxRetries    = 1
    errorStrategy = { task.exitStatus in 130..140 ? 'retry' : 'finish' }
  }
  withName: 'picard' {
    memory = { memForProcess(params.genomeSize, 1.0, 2, task.attempt) }
  }
  withName: 'gatk' {
    memory = { memForProcess(params.genomeSize, 1.0, 2, task.attempt) }
  }
}

profiles {
```

- [ ] **Step 2: Verify Nextflow parses the config without errors**

Run from the repo root:
```bash
nextflow config -profile processSingleExperiment
```
Expected: config dump with no errors. Confirm `genomeSize = 0` appears in params output.

- [ ] **Step 3: Commit**

```bash
git add nextflow.config
git commit -m "Add genome-size-based memory closures to nextflow.config"
```

---

### Task 2: Remove hardcoded `memory` directives from `alignment.nf`

**Files:**
- Modify: `modules/alignment.nf:107` (picard) and `modules/alignment.nf:147` (gatk)

The `jvmMem` derivation (`task.memory.toGiga() - 1`) in both process scripts must stay — it reads from `task.memory` which is now populated by the config closure. Only the `memory { ... }` directive line itself is removed.

- [ ] **Step 1: Remove `memory` directive from `picard` process**

In `modules/alignment.nf`, find the `picard` process block and remove this line:
```groovy
  memory { 5.GB * task.attempt }
```

The block should go from:
```groovy
process picard {
  container 'broadinstitute/picard:2.25.0'
  memory { 5.GB * task.attempt }

  input:
```
to:
```groovy
process picard {
  container 'broadinstitute/picard:2.25.0'

  input:
```

- [ ] **Step 2: Remove `memory` directive from `gatk` process**

In `modules/alignment.nf`, find the `gatk` process block and remove this line:
```groovy
  memory { 5.GB * task.attempt }
```

The block should go from:
```groovy
process gatk {
  container 'broadinstitute/gatk3:3.8-1'
  memory { 5.GB * task.attempt }

  input:
```
to:
```groovy
process gatk {
  container 'broadinstitute/gatk3:3.8-1'

  input:
```

- [ ] **Step 3: Verify pipeline still runs with stub tasks**

```bash
nextflow run main.nf -entry processSingleExperiment -profile processSingleExperiment -stub
```
Expected: all processes show as STUB COMPLETED, no config errors.

- [ ] **Step 4: Commit**

```bash
git add modules/alignment.nf
git commit -m "Remove hardcoded memory directives from picard and gatk"
```

---

### Task 3: Simplify `MakeDnaSeqNextflowConfig.pm`

**Files:**
- Modify: `ApiCommonWorkflow/Main/lib/perl/WorkflowSteps/MakeDnaSeqNextflowConfig.pm`

Keep the genome size calculation loop (lines 49–63). Remove `$bwaDefaultMemMb`, `$bwaRetryMemMb`, `$bwaRetries`, `$bwaBlock`, `$freebayesBlock`, and `$maxMemoryGigs` entirely. Emit `genomeSize = <bytes>` in the params block. Simplify the generated `process` block to just `executor` and `queue`.

- [ ] **Step 1: Remove `$maxMemoryGigs` fetch and `$bwaBlock`/`$freebayesBlock` variables**

Replace the section from line 41 to line 98 (everything after `my $bwaThreads` through `$bwaBlock .= "  }\n";`) with just the genome size calculation and a single `$genomeSizeBytes` variable:

```perl
    my $executor        = $self->getClusterExecutor();
    my $queue           = $self->getClusterQueue();
    my $genomeFastaFile = $self->getWorkflowDataDir() . "/" . $self->getParamValue("genomeFastaFile");

    my $genomeSizeBytes = 0;
    if (defined($genomeFastaFile) && -e $genomeFastaFile) {
        open(my $fh, "<", $genomeFastaFile) or die "Cannot open genome fasta '$genomeFastaFile': $!";
        while (<$fh>) {
            next if /^>/;
            chomp;
            $genomeSizeBytes += length($_);
        }
        close($fh);
    }
```

- [ ] **Step 2: Update the generated config printed in the `run` method**

Replace the `print F` block (lines 103–130) with:

```perl
        open(F, ">", $nextflowConfigFile) or die "$! :Can't open config file '$nextflowConfigFile' for writing";
        print F "
params {
  samplesheet              = \"$digestedSampleSheet\"
  bwaThreads               = $bwaThreads
  minCoverage              = $minCoverage
  genomeFastaFile          = \"$digestedGenomeFile\"
  gtfFile                  = \"$digestedGtfFile\"
  footprintFile            = \"$digestedFootprintFile\"
  winLen                   = $winLen
  ploidy                   = $ploidy
  outputDir                = \"$digestedResultsDir\"
  geneSourceIdOrthologFile = \"$digestedOrthologFile\"
  chrsForCalcFile          = \"$digestedChrsForCalcFile\"
  genomeSize               = $genomeSizeBytes
}

process {
  executor = '$executor'
  queue    = '$queue'
}

singularity {
  enabled    = true
  autoMounts = true
}
";
        close(F);
```

- [ ] **Step 3: Verify the final state of the file looks correct**

The `run` sub should now have:
1. All `getParamValue` calls (samplesheet, genomeFile, etc.) — unchanged
2. All `relativePathToNextflowClusterPath` calls — unchanged
3. `getConfig` calls for `minCoverage`, `winLen`, `bwaThreads` — unchanged
4. `getClusterExecutor` and `getClusterQueue` — unchanged
5. Genome FASTA loop producing `$genomeSizeBytes` — kept
6. `$isLsf`, `$bwaBlock`, `$freebayesBlock`, `$maxMemoryGigs`, `$bwaDefaultMemMb`, `$bwaRetryMemMb`, `$bwaRetries` — all gone
7. Generated config: params block includes `genomeSize`, process block has only `executor` and `queue`

Confirm by reading the file top-to-bottom and checking no dead variables remain.

- [ ] **Step 4: Commit**

```bash
cd /home/jbrestel/workspaces/dataLoad/project_home/ApiCommonWorkflow
git add Main/lib/perl/WorkflowSteps/MakeDnaSeqNextflowConfig.pm
git commit -m "Pass genomeSize param; remove memory math from MakeDnaSeqNextflowConfig"
```

---

### Task 4: End-to-end verification

- [ ] **Step 1: Confirm `nextflow config` shows expected memory for a known genome size**

```bash
nextflow config -profile processSingleExperiment
```

With `genomeSize = 0` (default), expected memory for all processes = `4 GB` (floor). Confirm the four `withName` blocks appear in output.

- [ ] **Step 2: Spot-check the formula manually for the test genome**

The test genome is at `data/genome_haploid/genome.fasta`. Get its size:
```bash
grep -v "^>" data/genome_haploid/genome.fasta | tr -d '\n' | wc -c
```

Then verify by hand:
- `genomeGb = <bytes> / 1_000_000_000`
- bwaMem: `round_up_to_pow2((genomeGb * 3.3 + 2) * 1.25)`
- FreeBayes: `round_up_to_pow2((genomeGb * 2.0 + 2) * 1.25)`
- Picard/GATK: `round_up_to_pow2((genomeGb * 1.0 + 2) * 1.25)`

For the Plasmodium falciparum 3D7 genome (~23 MB), all four processes should floor to **4 GB**.

- [ ] **Step 3: Run the stub pipeline to confirm no config errors**

```bash
nextflow run main.nf -entry processSingleExperiment -profile processSingleExperiment -stub
```
Expected: all processes STUB COMPLETED.

- [ ] **Step 4: Run the test suite**

```bash
nextflow run main.nf -entry runTests -profile tests
```
Expected: all Perl tests pass.
