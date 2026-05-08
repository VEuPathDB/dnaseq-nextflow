# VCF + Coverage BED Consensus FASTA Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace FreeBayes gVCF output with an explicit per-sample coverage BED, and rebuild consensus FASTA generation from VCF + BED instead of gVCF span rows.

**Architecture:** FreeBayes drops `--gvcf` and emits variant-only VCF. A new `makeCoverageBed` process derives covered regions (depth >= minCoverage) from the BAM via `bedtools genomecov`. A new Python script builds consensus FASTA by walking sparse VCF records and filling gaps using the coverage BED. A new `mergeCoverageBeds` process publishes a multi-sample coverage TSV using `bedtools multiinter`.

**Tech Stack:** Python 3 + cyvcf2, Nextflow DSL2, bedtools, bgzip/tabix, pytest

---

## File Map

| File | Action |
|---|---|
| `bin/makeConsensusFastaFromVcfAndBed.py` | **Create** — new consensus FASTA script |
| `bin/makeConsensusFastaFromGvcf.py` | **Delete** |
| `bin/splitGvcfAtZeroCoverage.py` | **Delete** |
| `testing/t/test_makeConsensusFastaFromVcfAndBed.py` | **Create** — pytest tests for new script |
| `modules/snp.nf` | **Modify** — remove gvcf from freebayes; remove `splitGvcfAtZeroCoverage` + `makeConsensusFromGvcf` + `mergeGvcfs` processes; add `makeCoverageBed` + `makeConsensusFromVcfAndBed` + `mergeCoverageBeds` |
| `modules/cnv.nf` | **Modify** — update `makeSnpDensity` input tuple (9 → 7 elements) |
| `workflows/processSingleExperiment.nf` | **Modify** — update imports, all freebayes tuple destructuring (9 → 7 elements), replace old process calls with new ones |

---

## Task 1: Write failing tests for `makeConsensusFastaFromVcfAndBed.py`

**Files:**
- Create: `testing/t/test_makeConsensusFastaFromVcfAndBed.py`

- [ ] **Step 1: Create the test file**

```python
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../bin'))

# This import will fail until the script is created — that's expected
from makeConsensusFastaFromVcfAndBed import (
    load_coverage_bed,
    fill_gap,
    is_covered,
    build_consensus,
    write_fasta,
)


class FakeVcfRecord:
    def __init__(self, chrom, pos, ref, alts, gt_bases_str):
        self.CHROM = chrom
        self.POS = pos        # 1-based
        self.REF = ref
        self.ALT = alts       # list of strings
        self.gt_bases = [gt_bases_str]  # single sample, e.g. 'A/T'
        self.FORMAT = 'GT'


class FakeVcf:
    def __init__(self, records_by_chrom):
        self._records = records_by_chrom

    def __call__(self, chrom):
        return iter(self._records.get(chrom, []))


# ── load_coverage_bed ─────────────────────────────────────────────────────────

def test_load_coverage_bed(tmp_path):
    bed = tmp_path / "sample.bed"
    bed.write_text("chr1\t0\t100\nchr1\t200\t300\nchr2\t50\t150\n")
    result = load_coverage_bed(str(bed))
    assert result == {'chr1': [(0, 100), (200, 300)], 'chr2': [(50, 150)]}


def test_load_coverage_bed_empty(tmp_path):
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    assert load_coverage_bed(str(bed)) == {}


# ── fill_gap ──────────────────────────────────────────────────────────────────

def test_fill_gap_all_covered():
    ref_seq = 'ACGTACGT'
    assert fill_gap(ref_seq, 0, 8, [(0, 8)]) == 'ACGTACGT'


def test_fill_gap_all_uncovered():
    assert fill_gap('ACGTACGT', 0, 8, []) == 'NNNNNNNN'


def test_fill_gap_partial_coverage():
    # covered [2,5), uncovered [0,2) and [5,8)
    ref_seq = 'ACGTACGT'
    result = fill_gap(ref_seq, 0, 8, [(2, 5)])
    assert result == 'NN' + 'GTA' + 'NNN'


def test_fill_gap_subrange_with_partial_overlap():
    # fill [1,6), coverage at [2,5)
    ref_seq = 'ACGTACGT'
    result = fill_gap(ref_seq, 1, 6, [(2, 5)])
    assert result == 'N' + 'GTA' + 'N'


def test_fill_gap_empty_range():
    assert fill_gap('ACGT', 3, 3, [(0, 4)]) == ''


def test_fill_gap_multiple_intervals():
    # covered [0,3) and [5,8), uncovered [3,5)
    ref_seq = 'ACGTACGT'
    result = fill_gap(ref_seq, 0, 8, [(0, 3), (5, 8)])
    assert result == 'ACG' + 'NN' + 'CGT'


# ── is_covered ────────────────────────────────────────────────────────────────

def test_is_covered_inside():
    assert is_covered([(0, 100)], 50) is True


def test_is_covered_at_start():
    assert is_covered([(0, 100)], 0) is True


def test_is_covered_at_end_exclusive():
    assert is_covered([(0, 100)], 100) is False


def test_is_covered_before_interval():
    assert is_covered([(50, 100)], 10) is False


def test_is_covered_empty():
    assert is_covered([], 0) is False


# ── build_consensus ───────────────────────────────────────────────────────────

def test_build_consensus_all_covered_no_variants():
    ref_seq = 'ACGTACGT'
    vcf = FakeVcf({})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(0, 8)])
    assert seq == 'ACGTACGT'


def test_build_consensus_all_uncovered_no_variants():
    ref_seq = 'ACGTACGT'
    vcf = FakeVcf({})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [])
    assert seq == 'NNNNNNNN'


def test_build_consensus_snp_covered():
    # SNP at pos 5 (1-based), A->T hom, fully covered
    ref_seq = 'ACGTACGT'
    record = FakeVcfRecord('chr1', 5, 'A', ['T'], 'T/T')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(0, 8)])
    assert seq == 'ACGT' + 'T' + 'CGT'


def test_build_consensus_snp_het_iupac():
    # het A/T at pos 1 → IUPAC W
    ref_seq = 'ACGT'
    record = FakeVcfRecord('chr1', 1, 'A', ['T'], 'A/T')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 4, ref_seq, vcf, [(0, 4)])
    assert seq[0] == 'W'
    assert seq[1:] == 'CGT'


def test_build_consensus_snp_uncovered_gives_n():
    ref_seq = 'ACGT'
    record = FakeVcfRecord('chr1', 1, 'A', ['T'], 'T/T')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 4, ref_seq, vcf, [])
    assert seq == 'NNNN'


def test_build_consensus_hom_insertion_covered():
    # Hom insertion at pos 3: REF=G, ALT=GCC (len 3 vs 1)
    ref_seq = 'ACGTACGT'
    record = FakeVcfRecord('chr1', 3, 'G', ['GCC'], 'GCC/GCC')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(0, 8)])
    # pos 0-1: AC, pos 2: G (covered gap = ref), variant emits GCC, then pos 3-7: TACGT
    assert seq == 'AC' + 'GCC' + 'TACGT'


def test_build_consensus_hom_deletion_covered():
    # Hom deletion at pos 1: REF=ACG, ALT=A (del of CG)
    ref_seq = 'ACGTACGT'
    record = FakeVcfRecord('chr1', 1, 'ACG', ['A'], 'A/A')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(0, 8)])
    # variant emits A (single base), then rest is TACGT
    assert seq == 'A' + 'TACGT'


def test_build_consensus_het_indel_01_emits_ref():
    # Het 0/1: one allele is REF → emit REF
    ref_seq = 'ACGTACGT'
    record = FakeVcfRecord('chr1', 2, 'CG', ['C'], 'CG/C')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(0, 8)])
    # pos 0: A, variant at pos 1 (0-based): REF=CG emitted, rest: TACGT
    assert seq == 'A' + 'CG' + 'TACGT'


def test_build_consensus_het_indel_12_emits_x_times_ref_len():
    # Het 1/2: both non-ref, REF=ACG (len 3) → XXX
    ref_seq = 'ACGTACGT'
    record = FakeVcfRecord('chr1', 1, 'ACG', ['AC', 'A'], 'AC/A')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(0, 8)])
    assert seq == 'XXX' + 'TACGT'


def test_build_consensus_gap_before_variant_partially_covered():
    # variant at pos 5 (1-based=0-based 4), coverage only [2,8)
    # gap [0,4): [0,2) uncovered → NN, [2,4) covered → ref GT
    ref_seq = 'ACGTACGT'
    record = FakeVcfRecord('chr1', 5, 'A', ['T'], 'T/T')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(2, 8)])
    assert seq == 'NN' + 'GT' + 'T' + 'CGT'


def test_build_consensus_gap_after_last_variant():
    # variant at pos 1, coverage [0,8), fill [1,8) after it
    ref_seq = 'ACGTACGT'
    record = FakeVcfRecord('chr1', 1, 'A', ['T'], 'T/T')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(0, 8)])
    assert seq == 'T' + 'CGTACGT'


def test_build_consensus_dot_gt_gives_n():
    ref_seq = 'ACGT'
    record = FakeVcfRecord('chr1', 1, 'A', ['T'], './.')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 4, ref_seq, vcf, [(0, 4)])
    assert seq[0] == 'N'


# ── write_fasta ───────────────────────────────────────────────────────────────

def test_write_fasta(tmp_path):
    out = tmp_path / "out.fa"
    with open(out, 'w') as fh:
        write_fasta(fh, 'chr1', 'ACGT' * 20)
    lines = out.read_text().splitlines()
    assert lines[0] == '>chr1'
    assert ''.join(lines[1:]) == 'ACGT' * 20
```

- [ ] **Step 2: Run tests to verify they all fail with ImportError**

```bash
cd /home/jbrestel/workspaces/dataLoad/project_home/dnaseq-nextflow
python -m pytest testing/t/test_makeConsensusFastaFromVcfAndBed.py -v 2>&1 | head -20
```

Expected: `ModuleNotFoundError: No module named 'makeConsensusFastaFromVcfAndBed'`

---

## Task 2: Implement `makeConsensusFastaFromVcfAndBed.py`

**Files:**
- Create: `bin/makeConsensusFastaFromVcfAndBed.py`

- [ ] **Step 1: Create the script**

```python
#!/usr/bin/env python3
import argparse
import bisect
import subprocess
from cyvcf2 import VCF

IUPAC = {
    frozenset('A'):    'A',
    frozenset('C'):    'C',
    frozenset('G'):    'G',
    frozenset('T'):    'T',
    frozenset('AC'):   'M',
    frozenset('AG'):   'R',
    frozenset('AT'):   'W',
    frozenset('CG'):   'S',
    frozenset('CT'):   'Y',
    frozenset('GT'):   'K',
    frozenset('ACG'):  'V',
    frozenset('ACT'):  'H',
    frozenset('AGT'):  'D',
    frozenset('CGT'):  'B',
    frozenset('ACGT'): 'N',
}


def get_chrom_seq(ref_fasta, chrom_name):
    result = subprocess.run(
        ['samtools', 'faidx', ref_fasta, chrom_name],
        capture_output=True, text=True, check=True,
    )
    return ''.join(result.stdout.splitlines()[1:]).upper()


def load_coverage_bed(bed_file):
    """Return dict[chrom] -> sorted list of (start, end) 0-based half-open intervals."""
    intervals = {}
    with open(bed_file) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            intervals.setdefault(chrom, []).append((start, end))
    return intervals


def fill_gap(ref_seq, gap_start, gap_end, chrom_intervals):
    """Fill [gap_start, gap_end) (0-based) with ref bases where covered, N elsewhere."""
    if gap_start >= gap_end:
        return ''
    if not chrom_intervals:
        return 'N' * (gap_end - gap_start)

    segments = []
    pos = gap_start
    starts = [iv[0] for iv in chrom_intervals]
    idx = max(bisect.bisect_right(starts, gap_start) - 1, 0)

    for iv_start, iv_end in chrom_intervals[idx:]:
        if iv_start >= gap_end:
            break
        fill_to = min(iv_start, gap_end)
        if pos < fill_to:
            segments.append('N' * (fill_to - pos))
            pos = fill_to
        cov_start = max(iv_start, pos)
        cov_end = min(iv_end, gap_end)
        if cov_start < cov_end:
            segments.append(ref_seq[cov_start:cov_end])
            pos = cov_end

    if pos < gap_end:
        segments.append('N' * (gap_end - pos))
    return ''.join(segments)


def is_covered(chrom_intervals, pos):
    """Return True if 0-based pos falls in any covered interval."""
    if not chrom_intervals:
        return False
    starts = [iv[0] for iv in chrom_intervals]
    idx = bisect.bisect_right(starts, pos) - 1
    if idx < 0:
        return False
    _, iv_end = chrom_intervals[idx]
    return pos < iv_end


def build_consensus(chrom_name, chrom_len, ref_seq, vcf, chrom_intervals):
    """
    Walk sparse VCF variant records for one chromosome and assemble consensus.

    Covered positions absent from VCF → reference bases (FreeBayes VCF-mode
    only emits variant sites; silence means ref).
    Positions outside coverage BED → N.
    SNPs        → IUPAC code.
    Hom indels  → actual allele sequence.
    Het 0/1     → REF (one allele is reference, no coordinate shift).
    Het 1/2     → X × len(REF) (ambiguous, both non-ref; indels.tsv excludes these).
    """
    segments = []
    ref_pos = 0  # 0-based

    for v in vcf(chrom_name):
        pos = v.POS - 1  # 1-based → 0-based

        if pos > ref_pos:
            segments.append(fill_gap(ref_seq, ref_pos, pos, chrom_intervals))
            ref_pos = pos

        ref_end = pos + len(v.REF)

        if not is_covered(chrom_intervals, pos):
            segments.append('N' * len(v.REF))
            ref_pos = ref_end
            continue

        if 'GT' not in v.FORMAT:
            segments.append('N' * len(v.REF))
            ref_pos = ref_end
            continue

        gt_str = v.gt_bases[0]
        if '.' in gt_str:
            segments.append('N' * len(v.REF))
            ref_pos = ref_end
            continue

        alleles = list(dict.fromkeys(gt_str.replace('|', '/').split('/')))

        if all(len(a) == 1 for a in alleles):
            segments.append(IUPAC.get(frozenset(alleles), 'N'))
            ref_pos = pos + 1
        elif len(alleles) == 1:
            segments.append(alleles[0])
            ref_pos = ref_end
        else:
            if v.REF in alleles:
                segments.append(v.REF)
            else:
                segments.append('X' * len(v.REF))
            ref_pos = ref_end

    if ref_pos < chrom_len:
        segments.append(fill_gap(ref_seq, ref_pos, chrom_len, chrom_intervals))

    return ''.join(segments)


def write_fasta(out_fh, name, seq, line_len=60):
    out_fh.write(f'>{name}\n')
    for i in range(0, len(seq), line_len):
        out_fh.write(seq[i:i + line_len] + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Build a consensus FASTA from a per-sample VCF and coverage BED.')
    parser.add_argument('-v', '--vcf',  required=True,
                        help='vcf.gz with tabix index (.tbi)')
    parser.add_argument('-b', '--bed',  required=True,
                        help='Coverage BED (regions with depth >= minCoverage)')
    parser.add_argument('-r', '--ref',  required=True,
                        help='Reference FASTA (indexed with samtools faidx)')
    parser.add_argument('-f', '--fai',  required=True,
                        help='Reference FASTA .fai index')
    parser.add_argument('-o', '--output', required=True)
    args = parser.parse_args()

    chroms = []
    with open(args.fai) as fh:
        for line in fh:
            parts = line.split('\t')
            chroms.append((parts[0], int(parts[1])))

    vcf = VCF(args.vcf)
    coverage = load_coverage_bed(args.bed)

    with open(args.output, 'w') as out:
        for chrom_name, chrom_len in chroms:
            ref_seq = get_chrom_seq(args.ref, chrom_name)
            chrom_intervals = coverage.get(chrom_name, [])
            seq = build_consensus(chrom_name, chrom_len, ref_seq, vcf, chrom_intervals)
            write_fasta(out, chrom_name, seq)


if __name__ == '__main__':
    main()
```

- [ ] **Step 2: Make executable**

```bash
chmod +x bin/makeConsensusFastaFromVcfAndBed.py
```

- [ ] **Step 3: Run tests and verify they pass**

```bash
cd /home/jbrestel/workspaces/dataLoad/project_home/dnaseq-nextflow
python -m pytest testing/t/test_makeConsensusFastaFromVcfAndBed.py -v
```

Expected: all tests PASS

- [ ] **Step 4: Commit**

```bash
git add bin/makeConsensusFastaFromVcfAndBed.py testing/t/test_makeConsensusFastaFromVcfAndBed.py
git commit -m "feat: add makeConsensusFastaFromVcfAndBed.py replacing gVCF-based consensus"
```

---

## Task 3: Update `freebayes` process — remove gVCF

**Files:**
- Modify: `modules/snp.nf:4-66`

- [ ] **Step 1: Remove `--gvcf` flag and gVCF output from the freebayes process**

In `modules/snp.nf`, replace the entire `freebayes` process with:

```groovy
process freebayes {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple val(sampleName), path(resultSortedGatkBam), path(resultSortedGatkBamIndex)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    tuple val(sampleName), path("${sampleName}.vcf.gz"), path("${sampleName}.vcf.gz.tbi"), path('freebayes.snps.vcf.gz'), path('freebayes.snps.vcf.gz.tbi'), path('freebayes.indels.vcf.gz'), path('freebayes.indels.vcf.gz.tbi'), emit: vcf_files

  script:
    """
    set -euo pipefail

    minAltFraction=\$([ "$params.ploidy" -eq 1 ] && echo "0.8" || echo "0.3")
    freebayes \\
      -f $genomeReorderedFasta \\
      -p $params.ploidy \\
      --min-coverage $params.minCoverage \\
      --min-alternate-fraction \$minAltFraction \\
      $resultSortedGatkBam | bcftools sort > freebayes.vcf

    bcftools view -v snps freebayes.vcf > freebayes.snps.vcf
    bcftools norm -m- freebayes.vcf | \\
      bcftools view --include 'strlen(ALT)!=strlen(REF) && ALT!~"^<"' > freebayes.indels.vcf

    bgzip freebayes.snps.vcf
    tabix -fp vcf freebayes.snps.vcf.gz
    bgzip freebayes.indels.vcf
    tabix -fp vcf freebayes.indels.vcf.gz

    bgzip freebayes.vcf
    mv freebayes.vcf.gz ${sampleName}.vcf.gz
    tabix -fp vcf ${sampleName}.vcf.gz
    """

  stub:
    """
    touch ${sampleName}.vcf.gz
    touch ${sampleName}.vcf.gz.tbi
    touch freebayes.snps.vcf.gz
    touch freebayes.snps.vcf.gz.tbi
    touch freebayes.indels.vcf.gz
    touch freebayes.indels.vcf.gz.tbi
    """
}
```

- [ ] **Step 2: Commit**

```bash
git add modules/snp.nf
git commit -m "feat: remove --gvcf from freebayes, drop gVCF output"
```

---

## Task 4: Add `makeCoverageBed`, `makeConsensusFromVcfAndBed`, and `mergeCoverageBeds` processes; remove old processes

**Files:**
- Modify: `modules/snp.nf`

- [ ] **Step 1: Remove `splitGvcfAtZeroCoverage`, `makeConsensusFromGvcf`, and `mergeGvcfs` process definitions**

Delete these three process blocks entirely from `modules/snp.nf` (lines 182-284 in the original — the `mergeGvcfs`, `splitGvcfAtZeroCoverage`, and `makeConsensusFromGvcf` blocks).

- [ ] **Step 2: Add the three new processes at the end of `modules/snp.nf`**

```groovy
process makeCoverageBed {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple val(sampleName), path(bamFile), path(bamIndex)

  output:
    tuple val(sampleName), path("${sampleName}.coverage.bed")

  script:
    """
    set -euo pipefail
    bedtools genomecov -ibam $bamFile -bga \\
      | awk -v mc=$params.minCoverage '\$4 >= mc' \\
      | cut -f1-3 \\
      | bedtools merge \\
      > ${sampleName}.coverage.bed
    """

  stub:
    """
    touch ${sampleName}.coverage.bed
    """
}


process makeConsensusFromVcfAndBed {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir", mode: "copy", saveAs: { "${sampleName}_consensus.fa.gz" }

  input:
    tuple val(sampleName), path(vcfGz), path(vcfGzTbi), path(coverageBed)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    path 'consensus.fa.gz'

  script:
    """
    set -euo pipefail
    makeConsensusFastaFromVcfAndBed.py \\
      --vcf $vcfGz \\
      --bed $coverageBed \\
      --ref $genomeReorderedFasta \\
      --fai $genomeReorderedFastaIndex \\
      --output consensus.fa
    bgzip consensus.fa
    """

  stub:
    """
    touch consensus.fa.gz
    """
}


process mergeCoverageBeds {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir", mode: "copy"

  input:
    path "*.coverage.bed"

  output:
    path 'coverage.tsv'

  script:
    """
    set -euo pipefail
    names=\$(ls *.coverage.bed | sed 's/\\.coverage\\.bed//' | tr '\\n' ' ')
    bedtools multiinter -names \$names -i *.coverage.bed > coverage.tsv
    """

  stub:
    """
    touch coverage.tsv
    """
}
```

- [ ] **Step 3: Commit**

```bash
git add modules/snp.nf
git commit -m "feat: add makeCoverageBed, makeConsensusFromVcfAndBed, mergeCoverageBeds processes"
```

---

## Task 5: Update `makeSnpDensity` input tuple in `cnv.nf`

**Files:**
- Modify: `modules/cnv.nf`

- [ ] **Step 1: Update the `makeSnpDensity` input to 7-element tuple (remove gvcf pair)**

Find the `makeSnpDensity` process input line:

```groovy
    tuple val(sampleName), path(freebayesVcfGz), path(freebayesVcfGzTbi), path(snpsVcfGz), path(snpsVcfGzTbi), path(indelsVcfGz), path(indelsVcfGzTbi), path(gvcfGz), path(gvcfGzTbi)
```

Replace with:

```groovy
    tuple val(sampleName), path(freebayesVcfGz), path(freebayesVcfGzTbi), path(snpsVcfGz), path(snpsVcfGzTbi), path(indelsVcfGz), path(indelsVcfGzTbi)
```

- [ ] **Step 2: Commit**

```bash
git add modules/cnv.nf
git commit -m "fix: update makeSnpDensity input tuple to 7 elements (no gvcf)"
```

---

## Task 6: Update `processSingleExperiment.nf` — imports, tuple maps, wiring

**Files:**
- Modify: `workflows/processSingleExperiment.nf`

- [ ] **Step 1: Update imports**

Replace these three import lines:

```groovy
include { splitGvcfAtZeroCoverage } from '../modules/snp.nf'
include { makeConsensusFromGvcf } from '../modules/snp.nf'
include { mergeGvcfs } from '../modules/snp.nf'
```

With:

```groovy
include { makeCoverageBed } from '../modules/snp.nf'
include { makeConsensusFromVcfAndBed } from '../modules/snp.nf'
include { mergeCoverageBeds } from '../modules/snp.nf'
```

- [ ] **Step 2: Update `combinedVcf` map (9 → 7 elements)**

Replace:

```groovy
    combinedVcf = freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi, gvcfGz, gvcfGzTbi ->
        tuple(sampleName, vcfGz, vcfGzTbi)
    }
```

With:

```groovy
    combinedVcf = freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, vcfGz, vcfGzTbi)
    }
```

- [ ] **Step 3: Update `makeIndelTSV` map (9 → 7 elements)**

Replace:

```groovy
    makeIndelTSV(freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi, gvcfGz, gvcfGzTbi ->
        tuple(sampleName, indelsVcfGz)
    }).collectFile(name: 'indels.tsv', storeDir: params.outputDir)
```

With:

```groovy
    makeIndelTSV(freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, indelsVcfGz)
    }).collectFile(name: 'indels.tsv', storeDir: params.outputDir)
```

- [ ] **Step 4: Replace gVCF-based consensus block with new coverage BED wiring**

Remove these lines (~106–121):

```groovy
    freebayesGvcf = freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi, gvcfGz, gvcfGzTbi ->
        tuple(sampleName, gvcfGz, gvcfGzTbi)
    }

    // Join gVCF with BAM so we can compute the full-genome BedGraph once and
    // use it to both split reference blocks at zero-coverage boundaries and
    // recompute per-sub-block DP values.
    splitGvcfResults = splitGvcfAtZeroCoverage(freebayesGvcf.join(gatkResults.bamTuple), reorderFastaResults)

    makeConsensusFromGvcf(splitGvcfResults, reorderFastaResults)

    mergeGvcfs(
        splitGvcfResults.count(),
        splitGvcfResults.map { sampleName, gvcfGz, gvcfGzTbi -> tuple(gvcfGz, gvcfGzTbi, "key") }
            .groupTuple(by: 2, sort: { a, b -> a <=> b })
    )
```

Replace with:

```groovy
    coverageBedResults = makeCoverageBed(gatkResults.bamTuple)

    perSampleVcf = freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, vcfGz, vcfGzTbi)
    }
    makeConsensusFromVcfAndBed(perSampleVcf.join(coverageBedResults), reorderFastaResults)

    mergeCoverageBeds(coverageBedResults.map { sampleName, bed -> bed }.collect())
```

- [ ] **Step 5: Update `snpsVcf` map inside `if (params.ploidy != 1)` block (9 → 7 elements)**

Replace:

```groovy
      snpsVcf = freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi, gvcfGz, gvcfGzTbi ->
        tuple(sampleName, snpsVcfGz, snpsVcfGzTbi)
      }
```

With:

```groovy
      snpsVcf = freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, snpsVcfGz, snpsVcfGzTbi)
      }
```

- [ ] **Step 6: Commit**

```bash
git add workflows/processSingleExperiment.nf
git commit -m "feat: wire makeCoverageBed + makeConsensusFromVcfAndBed + mergeCoverageBeds in processSingleExperiment"
```

---

## Task 7: Delete obsolete scripts

**Files:**
- Delete: `bin/makeConsensusFastaFromGvcf.py`
- Delete: `bin/splitGvcfAtZeroCoverage.py`

- [ ] **Step 1: Delete the scripts**

```bash
git rm bin/makeConsensusFastaFromGvcf.py bin/splitGvcfAtZeroCoverage.py
```

- [ ] **Step 2: Commit**

```bash
git commit -m "chore: delete makeConsensusFastaFromGvcf.py and splitGvcfAtZeroCoverage.py"
```

---

## Task 8: Stub-run verification

- [ ] **Step 1: Run the pipeline in stub mode to verify wiring**

```bash
cd /home/jbrestel/workspaces/dataLoad/project_home/dnaseq-nextflow
nextflow run main.nf -profile processSingleExperiment -stub-run 2>&1 | tail -30
```

Expected: workflow completes without errors. All processes resolve. No "No such variable" or tuple-size mismatch errors.

- [ ] **Step 2: If stub-run passes, run the Python tests one final time**

```bash
python -m pytest testing/t/test_makeConsensusFastaFromVcfAndBed.py -v
```

Expected: all tests PASS
