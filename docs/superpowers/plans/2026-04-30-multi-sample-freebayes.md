# Multi-Sample FreeBayes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace per-sample FreeBayes → bcftools merge with joint multi-sample FreeBayes over genome chunks, eliminating merge artifacts and ensuring coverage-consistent consensus FASTAs.

**Architecture:** Split the reference into fixed-size windows (`makeRegionBed`), run `freebayesMultiSample` per region, compute per-region multi-sample zero-coverage BED, split gVCF chunks at those boundaries, then gather into a single joint multi-sample VCF from which consensus FASTAs and per-sample CNV VCFs are derived.

**Tech Stack:** Nextflow DSL2, FreeBayes, bcftools, bedtools, samtools, Python 3 (cyvcf2), Perl

---

## File Map

| Action | File | What changes |
|---|---|---|
| Modify | `nextflow.config` | Add `chunkSize` param |
| Modify | `modules/snp.nf` | Add `makeRegionBed`, `makeMultiSampleZeroCoverageBed`, `concatMultiSampleVcf`, `extractSampleVcf`; modify `freebayesMultiSample`, `splitGvcfAtZeroCoverage`, `makeConsensusFromGvcf`, `makeIndelTSV`; remove `freebayes`, `mergeVcfs`, `makeMergedVariantIndex`, `mergeGvcfs`, `bcftoolsMpileupGvcf` |
| Modify | `modules/cnv.nf` | Shrink `makeSnpDensity` input tuple |
| Modify | `bin/splitGvcfAtZeroCoverage.py` | Replace `--bedgraph` with `--zero-cov-bed`; remove DP recomputation |
| Modify | `bin/makeConsensusFastaFromGvcf.py` | Handle multi-sample gVCF; output one FASTA per sample |
| Modify | `workflows/processSingleExperiment.nf` | Rewire all SNP channel logic |
| Create | `testing/t/test_splitGvcfAtZeroCoverage.py` | Pytest unit tests |
| Create | `testing/t/test_makeConsensusFastaFromGvcf.py` | Pytest unit tests |

---

### Task 1: Add `chunkSize` param and `makeRegionBed` process

**Files:**
- Modify: `nextflow.config`
- Modify: `modules/snp.nf`

- [ ] **Step 1: Add `chunkSize` to `nextflow.config` inside the `processSingleExperiment` params block (after `winLen = 1000`)**

```groovy
chunkSize = 1000000
```

- [ ] **Step 2: Add `makeRegionBed` to `modules/snp.nf` after the `freebayesMultiSample` process**

```nextflow
process makeRegionBed {
  container 'veupathdb/shortreadaligner:1.0.0'

  input:
    tuple path(genomeFasta), path(genomeFastaIndex)

  output:
    path 'regions.bed'

  script:
    """
    set -euo pipefail
    bedtools makewindows -g $genomeFastaIndex -w $params.chunkSize > regions.bed
    """

  stub:
    """
    printf 'chr1\\t0\\t1000000\\n' > regions.bed
    printf 'chr1\\t1000000\\t2000000\\n' >> regions.bed
    """
}
```

- [ ] **Step 3: Verify the file parses with no Nextflow syntax errors**

```bash
cd /home/jbrestel/workspaces/dataLoad/project_home/dnaseq-nextflow
nextflow run main.nf -entry processSingleExperiment -profile processSingleExperiment -stub 2>&1 | head -30
```

Expected: warnings about stub outputs are fine; a `Caused by` parse error is not.

- [ ] **Step 4: Commit**

```bash
git add nextflow.config modules/snp.nf
git commit -m "feat: add makeRegionBed process and chunkSize param"
```

---

### Task 2: Add `makeMultiSampleZeroCoverageBed` process

**Files:**
- Modify: `modules/snp.nf`

- [ ] **Step 1: Add `makeMultiSampleZeroCoverageBed` to `modules/snp.nf` after `makeRegionBed`**

```nextflow
process makeMultiSampleZeroCoverageBed {
  container 'veupathdb/shortreadaligner:1.0.0'

  input:
    path bamFiles
    path bamIndexes
    val regionLine

  output:
    tuple val(regionKey), path('zero_cov.bed')

  script:
    def fields = regionLine.tokenize('\t')
    def chrom  = fields[0]
    def start  = fields[1].toLong()
    def end    = fields[2].toLong()
    regionKey  = "${chrom}_${start}_${end}"
    """
    set -euo pipefail
    for bam in *.bam; do
      samtools view -b -h \$bam ${chrom}:$((start + 1))-${end} > region_\${bam%.bam}.bam
      samtools index region_\${bam%.bam}.bam
      bedtools genomecov -ibam region_\${bam%.bam}.bam -bga | \\
        awk '\$4 == 0 {print \$1 "\\t" \$2 "\\t" \$3}' >> all_zero.bed
    done
    if [ -s all_zero.bed ]; then
      bedtools sort -i all_zero.bed | bedtools merge > zero_cov.bed
    else
      touch zero_cov.bed
    fi
    """

  stub:
    def fields = regionLine.tokenize('\t')
    def chrom  = fields[0]
    def start  = fields[1].toLong()
    def end    = fields[2].toLong()
    regionKey  = "${chrom}_${start}_${end}"
    """
    touch zero_cov.bed
    """
}
```

- [ ] **Step 2: Verify parse**

```bash
nextflow run main.nf -entry processSingleExperiment -profile processSingleExperiment -stub 2>&1 | grep -i "error\|exception" | head -10
```

Expected: no parse errors.

- [ ] **Step 3: Commit**

```bash
git add modules/snp.nf
git commit -m "feat: add makeMultiSampleZeroCoverageBed process"
```

---

### Task 3: Modify `freebayesMultiSample` to accept a region and emit a region key

**Files:**
- Modify: `modules/snp.nf`

- [ ] **Step 1: Replace the existing `freebayesMultiSample` process with this version**

```nextflow
process freebayesMultiSample {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    path bamFiles
    path bamIndexes
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)
    val regionLine

  output:
    tuple val(regionKey), path("${regionKey}.g.vcf.gz"), path("${regionKey}.g.vcf.gz.tbi"), emit: gvcf

  script:
    def fields    = regionLine.tokenize('\t')
    def chrom     = fields[0]
    def start     = fields[1].toLong()
    def end       = fields[2].toLong()
    regionKey     = "${chrom}_${start}_${end}"
    def regionStr = "${chrom}:${start}-${end}"
    """
    set -euo pipefail
    ls *.bam > bam_list.txt
    minAltFraction=\$([ "$params.ploidy" -eq 1 ] && echo "0.8" || echo "0.3")
    freebayes \\
      -f $genomeReorderedFasta \\
      -p $params.ploidy \\
      --min-coverage $params.minCoverage \\
      --min-alternate-fraction \$minAltFraction \\
      --gvcf \\
      --region ${regionStr} \\
      --bam-list bam_list.txt \\
    | bcftools sort -O z -o ${regionKey}.g.vcf.gz
    bcftools index -t ${regionKey}.g.vcf.gz
    """

  stub:
    def fields = regionLine.tokenize('\t')
    def chrom  = fields[0]
    def start  = fields[1].toLong()
    def end    = fields[2].toLong()
    regionKey  = "${chrom}_${start}_${end}"
    """
    touch ${regionKey}.g.vcf.gz ${regionKey}.g.vcf.gz.tbi
    """
}
```

- [ ] **Step 2: Verify parse**

```bash
nextflow run main.nf -entry processSingleExperiment -profile processSingleExperiment -stub 2>&1 | grep -i "error\|exception" | head -10
```

- [ ] **Step 3: Commit**

```bash
git add modules/snp.nf
git commit -m "feat: add region input and key output to freebayesMultiSample"
```

---

### Task 4: Write tests then modify `splitGvcfAtZeroCoverage.py`

**Files:**
- Create: `testing/t/test_splitGvcfAtZeroCoverage.py`
- Modify: `bin/splitGvcfAtZeroCoverage.py`

- [ ] **Step 1: Write the failing tests**

Create `testing/t/test_splitGvcfAtZeroCoverage.py`:

```python
import os
import sys
import tempfile
import textwrap

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../bin'))

from splitGvcfAtZeroCoverage import load_zero_cov_bed, covered_intervals, process_gvcf


def test_load_zero_cov_bed_parses_intervals():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write('chr1\t100\t200\n')
        f.write('chr1\t500\t600\n')
        f.write('chr2\t0\t50\n')
        fname = f.name
    try:
        result = load_zero_cov_bed(fname)
        assert result['chr1'] == [(100, 200), (500, 600)]
        assert result['chr2'] == [(0, 50)]
    finally:
        os.unlink(fname)


def test_load_zero_cov_bed_empty_file():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        fname = f.name
    try:
        result = load_zero_cov_bed(fname)
        assert result == {}
    finally:
        os.unlink(fname)


def test_process_gvcf_splits_ref_block_at_zero_coverage():
    """A reference block spanning a zero-coverage gap should be split into two sub-blocks."""
    gvcf_content = textwrap.dedent("""\
        ##fileformat=VCFv4.2
        ##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
        chr1\t1\t.\tA\t<*>\t.\t.\tEND=1000\tGT:DP\t0/0:20\t0/0:15
    """)
    bed_content = 'chr1\t499\t501\n'  # zero-coverage gap at positions 500-501

    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as gvcf_f, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as bed_f, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as out_f:
        gvcf_f.write(gvcf_content)
        bed_f.write(bed_content)
        gvcf_path = gvcf_f.name
        bed_path = bed_f.name
        out_path = out_f.name

    try:
        process_gvcf(gvcf_path, bed_path, out_path)
        with open(out_path) as f:
            lines = [l for l in f if not l.startswith('#')]
        assert len(lines) == 2, f"Expected 2 sub-blocks, got {len(lines)}: {lines}"
        pos1 = int(lines[0].split('\t')[1])
        pos2 = int(lines[1].split('\t')[1])
        assert pos1 == 1
        assert pos2 == 502  # first base after zero-cov gap (501 0-based = 502 1-based)
    finally:
        for p in [gvcf_path, bed_path, out_path]:
            os.unlink(p)


def test_process_gvcf_drops_fully_zero_block():
    """A reference block entirely within a zero-coverage region should be dropped."""
    gvcf_content = textwrap.dedent("""\
        ##fileformat=VCFv4.2
        ##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
        chr1\t100\t.\tA\t<*>\t.\t.\tEND=200\tGT:DP\t0/0:0
    """)
    bed_content = 'chr1\t99\t201\n'  # covers entire block (0-based 99-201 = 1-based 100-201)

    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as gvcf_f, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as bed_f, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as out_f:
        gvcf_f.write(gvcf_content)
        bed_f.write(bed_content)
        gvcf_path = gvcf_f.name
        bed_path = bed_f.name
        out_path = out_f.name

    try:
        process_gvcf(gvcf_path, bed_path, out_path)
        with open(out_path) as f:
            lines = [l for l in f if not l.startswith('#')]
        assert len(lines) == 0, f"Expected block to be dropped, got: {lines}"
    finally:
        for p in [gvcf_path, bed_path, out_path]:
            os.unlink(p)


def test_process_gvcf_passes_variant_records_unchanged():
    """Non-reference-block records should pass through untouched."""
    gvcf_content = textwrap.dedent("""\
        ##fileformat=VCFv4.2
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
        chr1\t100\t.\tA\tT\t50\t.\t.\tGT:DP\t1/1:30\t0/1:25
    """)
    bed_content = 'chr1\t99\t101\n'

    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as gvcf_f, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as bed_f, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as out_f:
        gvcf_f.write(gvcf_content)
        bed_f.write(bed_content)
        gvcf_path = gvcf_f.name
        bed_path = bed_f.name
        out_path = out_f.name

    try:
        process_gvcf(gvcf_path, bed_path, out_path)
        with open(out_path) as f:
            lines = [l for l in f if not l.startswith('#')]
        assert len(lines) == 1
        assert lines[0].split('\t')[1] == '100'
    finally:
        for p in [gvcf_path, bed_path, out_path]:
            os.unlink(p)
```

- [ ] **Step 2: Run the tests — verify they fail because the new API doesn't exist yet**

```bash
cd /home/jbrestel/workspaces/dataLoad/project_home/dnaseq-nextflow
python -m pytest testing/t/test_splitGvcfAtZeroCoverage.py -v 2>&1 | tail -20
```

Expected: `ImportError` or `AttributeError` on `load_zero_cov_bed` since it doesn't exist yet.

- [ ] **Step 3: Replace `bin/splitGvcfAtZeroCoverage.py` with the new version**

Key changes: replace `load_bedgraph` with `load_zero_cov_bed`; remove `compute_depth_stats` and `update_dp`; update `process_gvcf` and `main`.

```python
#!/usr/bin/env python3
"""
Split gVCF reference blocks (<*>) at zero-coverage boundaries.

Reads a pre-computed zero-coverage BED file (union of all-sample zero-coverage
intervals, produced by makeMultiSampleZeroCoverageBed) to identify where
reference blocks should be split or dropped.

Variant records are passed through unchanged.
Works with both single-sample and multi-sample gVCFs.
"""

import sys
import gzip
import argparse
import bisect
from collections import defaultdict


def load_fasta(fasta_file):
    """Load a FASTA file into a dict of chrom -> uppercase sequence string."""
    sequences = {}
    current_chrom = None
    chunks = []
    opener = gzip.open if fasta_file.endswith('.gz') else open
    with opener(fasta_file, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if current_chrom is not None:
                    sequences[current_chrom] = ''.join(chunks)
                current_chrom = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
    if current_chrom is not None:
        sequences[current_chrom] = ''.join(chunks)
    return sequences


def load_zero_cov_bed(bed_file):
    """
    Load a BED file of zero-coverage intervals.
    Returns: dict[chrom] -> list of (start, end)  (0-based half-open, sorted)
    """
    zero_regions = defaultdict(list)
    with open(bed_file) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            zero_regions[chrom].append((start, end))
    return zero_regions


def covered_intervals(pos1, end1, zero_regions_chrom):
    """
    Given a block [pos1, end1] (1-based inclusive) and zero-coverage intervals
    (0-based half-open, sorted), return list of covered (start, end) tuples
    (1-based inclusive).
    """
    block_start = pos1 - 1  # 0-based
    block_end = end1        # 0-based half-open

    covered = []
    current = block_start

    for z_start, z_end in zero_regions_chrom:
        if z_end <= block_start:
            continue
        if z_start >= block_end:
            break
        if current < z_start:
            covered.append((current + 1, z_start))  # back to 1-based
        current = max(current, z_end)

    if current < block_end:
        covered.append((current + 1, block_end))  # back to 1-based

    return covered


def update_end_in_info(info, new_end):
    fields = info.split(';')
    return ';'.join('END={}'.format(new_end) if f.startswith('END=') else f for f in fields)


def process_gvcf(gvcf_file, zero_cov_bed_file, output_file, ref_fasta=None):
    zero_regions = load_zero_cov_bed(zero_cov_bed_file)
    ref_seqs = load_fasta(ref_fasta) if ref_fasta else {}

    in_opener = gzip.open if gvcf_file.endswith('.gz') else open

    with in_opener(gvcf_file, 'rt') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue

            parts = line.rstrip('\n').split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            alt = parts[4]
            info = parts[7]

            # Pass variant records through unchanged
            if alt != '<*>':
                fout.write(line)
                continue

            # Parse END from INFO
            end = pos
            for field in info.split(';'):
                if field.startswith('END='):
                    end = int(field[4:])
                    break

            # Clamp END to chromosome length
            if chrom in ref_seqs:
                end = min(end, len(ref_seqs[chrom]))

            chrom_zeros = zero_regions.get(chrom, [])
            if not chrom_zeros:
                fout.write(line)
                continue

            intervals = covered_intervals(pos, end, chrom_zeros)

            if not intervals:
                # Entire block is zero-coverage; drop it
                continue

            if len(intervals) == 1 and intervals[0] == (pos, end):
                # No splitting needed
                fout.write(line)
                continue

            # Emit one record per covered sub-interval
            for sub_start, sub_end in intervals:
                new_parts = parts[:]
                new_parts[1] = str(sub_start)
                new_parts[7] = update_end_in_info(info, sub_end)
                if sub_start != pos and chrom in ref_seqs:
                    new_parts[3] = ref_seqs[chrom][sub_start - 1]
                fout.write('\t'.join(new_parts) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gvcf',         required=True,
                        help='Input gVCF (plain or bgzipped)')
    parser.add_argument('--zero-cov-bed', required=True, dest='zero_cov_bed',
                        help='BED file of zero-coverage intervals (union across all samples)')
    parser.add_argument('--output',       required=True,
                        help='Output gVCF')
    parser.add_argument('--ref',          required=False,
                        help='Reference FASTA; fixes REF alleles on split sub-blocks')
    args = parser.parse_args()
    process_gvcf(args.gvcf, args.zero_cov_bed, args.output, ref_fasta=args.ref)


if __name__ == '__main__':
    main()
```

- [ ] **Step 4: Run the tests — verify they all pass**

```bash
python -m pytest testing/t/test_splitGvcfAtZeroCoverage.py -v 2>&1 | tail -20
```

Expected: 4 tests PASSED.

- [ ] **Step 5: Commit**

```bash
git add bin/splitGvcfAtZeroCoverage.py testing/t/test_splitGvcfAtZeroCoverage.py
git commit -m "feat: replace --bedgraph with --zero-cov-bed in splitGvcfAtZeroCoverage; support multi-sample gVCF"
```

---

### Task 5: Modify the `splitGvcfAtZeroCoverage` Nextflow process

**Files:**
- Modify: `modules/snp.nf`

- [ ] **Step 1: Replace the `splitGvcfAtZeroCoverage` process**

```nextflow
process splitGvcfAtZeroCoverage {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple val(regionKey), path(gvcfGz, stageAs: 'input.g.vcf.gz'), path(gvcfGzTbi, stageAs: 'input.g.vcf.gz.tbi'), path(zeroCovBed)
    tuple path(genomeFasta), path(genomeFastaIndex)

  output:
    tuple val(regionKey), path("${regionKey}.split.g.vcf.gz"), path("${regionKey}.split.g.vcf.gz.tbi")

  script:
    """
    set -euo pipefail
    splitGvcfAtZeroCoverage.py \\
      --gvcf input.g.vcf.gz \\
      --zero-cov-bed $zeroCovBed \\
      --ref $genomeFasta \\
      --output /dev/stdout \\
    | bcftools sort -O z -o ${regionKey}.split.g.vcf.gz
    bcftools index -t ${regionKey}.split.g.vcf.gz
    """

  stub:
    """
    touch ${regionKey}.split.g.vcf.gz ${regionKey}.split.g.vcf.gz.tbi
    """
}
```

- [ ] **Step 2: Verify parse**

```bash
nextflow run main.nf -entry processSingleExperiment -profile processSingleExperiment -stub 2>&1 | grep -i "error\|exception" | head -10
```

- [ ] **Step 3: Commit**

```bash
git add modules/snp.nf
git commit -m "feat: update splitGvcfAtZeroCoverage process to use zero-cov-bed and emit region key"
```

---

### Task 6: Add `concatMultiSampleVcf` process

**Files:**
- Modify: `modules/snp.nf`

- [ ] **Step 1: Add `concatMultiSampleVcf` to `modules/snp.nf` after `splitGvcfAtZeroCoverage`**

```nextflow
process concatMultiSampleVcf {
  container 'veupathdb/shortreadaligner:1.0.0'

  publishDir "$params.outputDir", mode: "copy"

  input:
    path gvcfFiles
    path gvcfIndexes

  output:
    tuple path('multisample.g.vcf.gz'), path('multisample.g.vcf.gz.tbi'), emit: gvcf
    tuple path('multisample.vcf.gz'),   path('multisample.vcf.gz.tbi'),   emit: vcf

  script:
    """
    set -euo pipefail
    ls *.split.g.vcf.gz | sort -V > file_list.txt
    bcftools concat --naive-force -f file_list.txt | bcftools sort -O z -o multisample.g.vcf.gz
    bcftools index -t multisample.g.vcf.gz
    bcftools view -e 'ALT[0]="<*>"' multisample.g.vcf.gz | bcftools sort -O z -o multisample.vcf.gz
    bcftools index -t multisample.vcf.gz
    """

  stub:
    """
    touch multisample.g.vcf.gz multisample.g.vcf.gz.tbi
    touch multisample.vcf.gz   multisample.vcf.gz.tbi
    """
}
```

- [ ] **Step 2: Commit**

```bash
git add modules/snp.nf
git commit -m "feat: add concatMultiSampleVcf process"
```

---

### Task 7: Write tests then modify `makeConsensusFastaFromGvcf.py`

**Files:**
- Create: `testing/t/test_makeConsensusFastaFromGvcf.py`
- Modify: `bin/makeConsensusFastaFromGvcf.py`

- [ ] **Step 1: Write the failing tests**

Create `testing/t/test_makeConsensusFastaFromGvcf.py`:

```python
import os
import sys
import tempfile
import textwrap

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../bin'))

from makeConsensusFastaFromGvcf import build_consensus, write_fasta

IUPAC = {
    frozenset('A'): 'A', frozenset('C'): 'C', frozenset('G'): 'G', frozenset('T'): 'T',
    frozenset('AC'): 'M', frozenset('AG'): 'R', frozenset('AT'): 'W',
    frozenset('CG'): 'S', frozenset('CT'): 'Y', frozenset('GT'): 'K',
}


class FakeVcfRecord:
    def __init__(self, chrom, pos, ref, alts, gt_bases_list, dp_list, end=None):
        self.CHROM = chrom
        self.POS = pos
        self.REF = ref
        self.ALT = alts
        self.gt_bases = gt_bases_list
        self._dp = dp_list
        self._end = end
        self.FORMAT = 'GT:DP'
        self.INFO = {}
        if end:
            self.INFO['END'] = end

    def format(self, key):
        if key == 'DP':
            return [[d] for d in self._dp]
        return None

    def INFO_get(self, key, default=None):
        return self.INFO.get(key, default)


class FakeVcf:
    def __init__(self, records_by_chrom):
        self._records = records_by_chrom
        self.samples = ['SAMPLE1', 'SAMPLE2']

    def __call__(self, chrom):
        return iter(self._records.get(chrom, []))


def _make_record_ref_block(chrom, pos, end, dp_s1, dp_s2):
    r = FakeVcfRecord(chrom, pos, 'A', None, ['A/A', 'A/A'], [dp_s1, dp_s2], end=end)
    r.INFO['END'] = end
    r.ALT = []
    return r


def _make_snp_record(chrom, pos, ref, alt, gt_s1, gt_s2, dp_s1=30, dp_s2=30):
    return FakeVcfRecord(chrom, pos, ref, [alt], [gt_s1, gt_s2], [dp_s1, dp_s2])


def test_build_consensus_ref_block_sample0():
    """REF block with sufficient coverage → reference bases in output."""
    ref_seq = 'ACGT' * 25  # 100 bp
    record = _make_record_ref_block('chr1', 1, 100, dp_s1=20, dp_s2=5)

    class MinimalVcf:
        def __call__(self, chrom):
            return iter([record])

    seq = build_consensus('chr1', 100, ref_seq, MinimalVcf(), min_coverage=1, sample_idx=0)
    assert seq == ref_seq


def test_build_consensus_ref_block_zero_coverage_masked():
    """REF block with dp < min_coverage → Ns."""
    ref_seq = 'A' * 50
    record = _make_record_ref_block('chr1', 1, 50, dp_s1=0, dp_s2=10)

    class MinimalVcf:
        def __call__(self, chrom):
            return iter([record])

    seq = build_consensus('chr1', 50, ref_seq, MinimalVcf(), min_coverage=1, sample_idx=0)
    assert seq == 'N' * 50


def test_build_consensus_snp_uses_correct_sample():
    """SNP record: each sample index yields the correct base/IUPAC code."""
    ref_seq = 'A' * 10
    record = _make_snp_record('chr1', 5, 'A', 'T', gt_s1='A/T', gt_s2='T/T')

    class MinimalVcf:
        def __call__(self, chrom):
            return iter([record])

    # sample 0: A/T het → IUPAC W
    seq0 = build_consensus('chr1', 10, ref_seq, MinimalVcf(), min_coverage=1, sample_idx=0)
    assert seq0[4] == 'W', f"Expected W (het A/T) for sample 0, got {seq0[4]}"

    # sample 1: T/T hom → T
    seq1 = build_consensus('chr1', 10, ref_seq, MinimalVcf(), min_coverage=1, sample_idx=1)
    assert seq1[4] == 'T', f"Expected T (hom T/T) for sample 1, got {seq1[4]}"


def test_write_fasta_produces_correct_output():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        fname = f.name
    try:
        with open(fname, 'w') as fh:
            write_fasta(fh, 'chr1', 'ACGT' * 20)
        with open(fname) as fh:
            lines = fh.read().splitlines()
        assert lines[0] == '>chr1'
        assert ''.join(lines[1:]) == 'ACGT' * 20
    finally:
        os.unlink(fname)
```

- [ ] **Step 2: Run tests — verify they fail**

```bash
python -m pytest testing/t/test_makeConsensusFastaFromGvcf.py -v 2>&1 | tail -20
```

Expected: failures because `build_consensus` doesn't yet accept `sample_idx`.

- [ ] **Step 3: Modify `bin/makeConsensusFastaFromGvcf.py`**

Change `build_consensus` to accept `sample_idx`, update `main` to loop over samples and write one file per sample.

The complete replacement for `bin/makeConsensusFastaFromGvcf.py`:

```python
#!/usr/bin/env python3
import argparse
import os
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
    """Load a single chromosome sequence from an indexed FASTA."""
    result = subprocess.run(
        ['samtools', 'faidx', ref_fasta, chrom_name],
        capture_output=True, text=True, check=True,
    )
    return ''.join(result.stdout.splitlines()[1:]).upper()


def build_consensus(chrom_name, chrom_len, ref_seq, vcf, min_coverage, sample_idx):
    """
    Walk through g.vcf records for one chromosome and assemble the consensus
    for the sample at sample_idx.
    """
    segments = []
    ref_pos = 0

    for v in vcf(chrom_name):
        pos = v.POS - 1  # VCF 1-based → 0-based

        if pos > ref_pos:
            segments.append('N' * (pos - ref_pos))
            ref_pos = pos

        dp_arr = v.format('DP')
        dp = int(dp_arr[sample_idx][0]) if dp_arr is not None else 0

        # REF block
        if not v.ALT or all(a == '.' or a.startswith('<') for a in v.ALT):
            end = (v.INFO.get('END') or v.POS) - 1  # 0-based inclusive
            span = end - pos + 1
            if dp >= min_coverage:
                segments.append(ref_seq[pos:end + 1])
            else:
                segments.append('N' * span)
            ref_pos = end + 1
            continue

        # Variant record
        if dp < min_coverage:
            segments.append('N' * len(v.REF))
            ref_pos = pos + len(v.REF)
            continue

        if 'GT' not in v.FORMAT:
            segments.append('N' * len(v.REF))
            ref_pos = pos + len(v.REF)
            continue

        gt_str = v.gt_bases[sample_idx]
        if '.' in gt_str:
            segments.append('N' * len(v.REF))
            ref_pos = pos + len(v.REF)
            continue

        alleles = list(dict.fromkeys(gt_str.replace('|', '/').split('/')))

        if all(len(a) == 1 for a in alleles):
            base = IUPAC.get(frozenset(alleles), 'N')
            segments.append(base)
            ref_pos = pos + 1

        elif len(alleles) == 1:
            segments.append(alleles[0])
            ref_pos = pos + len(v.REF)

        else:
            if v.REF in alleles:
                segments.append(v.REF)
            else:
                segments.append('X' * len(v.REF))
            ref_pos = pos + len(v.REF)

    if ref_pos < chrom_len:
        segments.append('N' * (chrom_len - ref_pos))

    return ''.join(segments)


def write_fasta(out_fh, name, seq, line_len=60):
    out_fh.write(f'>{name}\n')
    for i in range(0, len(seq), line_len):
        out_fh.write(seq[i:i + line_len] + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Build per-sample consensus FASTAs from a multi-sample g.vcf.')
    parser.add_argument('-g', '--gvcf',          required=True,
                        help='g.vcf.gz with tabix index (.tbi)')
    parser.add_argument('-r', '--ref',           required=True,
                        help='Reference FASTA (indexed with samtools faidx)')
    parser.add_argument('-f', '--fai',           required=True,
                        help='Reference FASTA .fai index')
    parser.add_argument('-mc', '--min-coverage', required=True, type=int,
                        dest='min_coverage')
    parser.add_argument('-o', '--output-dir',    required=False, default='.',
                        dest='output_dir',
                        help='Directory for output FASTA files (default: current directory)')
    args = parser.parse_args()

    chroms = []
    with open(args.fai) as fh:
        for line in fh:
            parts = line.split('\t')
            chroms.append((parts[0], int(parts[1])))

    # Read sample names from header, then iterate once per sample
    header_vcf = VCF(args.gvcf)
    sample_names = list(header_vcf.samples)

    for sample_idx, sample_name in enumerate(sample_names):
        vcf = VCF(args.gvcf)
        out_path = os.path.join(args.output_dir, f'{sample_name}_consensus.fa')
        with open(out_path, 'w') as out:
            for chrom_name, chrom_len in chroms:
                ref_seq = get_chrom_seq(args.ref, chrom_name)
                seq = build_consensus(chrom_name, chrom_len, ref_seq, vcf,
                                      args.min_coverage, sample_idx)
                write_fasta(out, chrom_name, seq)


if __name__ == '__main__':
    main()
```

- [ ] **Step 4: Run tests — verify they all pass**

```bash
python -m pytest testing/t/test_makeConsensusFastaFromGvcf.py -v 2>&1 | tail -20
```

Expected: 4 tests PASSED.

- [ ] **Step 5: Commit**

```bash
git add bin/makeConsensusFastaFromGvcf.py testing/t/test_makeConsensusFastaFromGvcf.py
git commit -m "feat: makeConsensusFastaFromGvcf supports multi-sample gVCF, outputs one FASTA per sample"
```

---

### Task 8: Modify `makeConsensusFromGvcf` Nextflow process

**Files:**
- Modify: `modules/snp.nf`

- [ ] **Step 1: Replace `makeConsensusFromGvcf` in `modules/snp.nf`**

```nextflow
process makeConsensusFromGvcf {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir", mode: "copy"

  input:
    tuple path(gvcfGz), path(gvcfGzTbi)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    path '*_consensus.fa.gz'

  script:
    """
    set -euo pipefail
    makeConsensusFastaFromGvcf.py \\
      --gvcf $gvcfGz \\
      --ref $genomeReorderedFasta \\
      --fai $genomeReorderedFastaIndex \\
      --min-coverage $params.minCoverage \\
      --output-dir .
    for f in *_consensus.fa; do
      bgzip \$f
    done
    """

  stub:
    """
    touch stub_consensus.fa.gz
    """
}
```

- [ ] **Step 2: Commit**

```bash
git add modules/snp.nf
git commit -m "feat: update makeConsensusFromGvcf process for multi-sample gVCF"
```

---

### Task 9: Modify `makeIndelTSV` to handle multi-sample gVCF

**Files:**
- Modify: `modules/snp.nf`

- [ ] **Step 1: Replace `makeIndelTSV` in `modules/snp.nf`**

```nextflow
process makeIndelTSV {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple path(gvcfGz), path(gvcfGzTbi)

  output:
    path('output.tsv')

  script:
    """
    set -euo pipefail
    bcftools query -l $gvcfGz > samples.txt

    while read sample; do
      bcftools view -s \$sample --trim-alt-alleles $gvcfGz | \\
        bcftools norm -m- | \\
        bcftools view --include 'strlen(ALT)!=strlen(REF) && ALT!~"^<"' | \\
        bgzip > \${sample}_indels.vcf.gz
      findValues.pl -i \${sample}_indels.vcf.gz -s \$sample -o \${sample}.tsv
    done < samples.txt

    cat *.tsv > output.tsv
    """

  stub:
    """
    touch output.tsv
    """
}
```

- [ ] **Step 2: Commit**

```bash
git add modules/snp.nf
git commit -m "feat: update makeIndelTSV to iterate over all samples in multi-sample gVCF"
```

---

### Task 10: Add `extractSampleVcf` process and update `makeSnpDensity`

**Files:**
- Modify: `modules/snp.nf`
- Modify: `modules/cnv.nf`

- [ ] **Step 1: Add `extractSampleVcf` to `modules/snp.nf` (after `makeIndelTSV`)**

```nextflow
process extractSampleVcf {
  container 'veupathdb/shortreadaligner:1.0.0'

  input:
    tuple val(sampleName), path(vcfGz), path(vcfGzTbi)

  output:
    tuple val(sampleName), path("${sampleName}.vcf.gz"), path("${sampleName}.vcf.gz.tbi"), path("${sampleName}.snps.vcf.gz"), path("${sampleName}.snps.vcf.gz.tbi"), path("${sampleName}.indels.vcf.gz"), path("${sampleName}.indels.vcf.gz.tbi"), emit: vcf_files

  script:
    """
    set -euo pipefail
    bcftools view -s $sampleName --trim-alt-alleles $vcfGz -O z -o ${sampleName}.vcf.gz
    bcftools index -t ${sampleName}.vcf.gz

    bcftools view -v snps ${sampleName}.vcf.gz -O z -o ${sampleName}.snps.vcf.gz
    bcftools index -t ${sampleName}.snps.vcf.gz

    bcftools norm -m- ${sampleName}.vcf.gz | \\
      bcftools view --include 'strlen(ALT)!=strlen(REF) && ALT!~"^<"' -O z -o ${sampleName}.indels.vcf.gz
    bcftools index -t ${sampleName}.indels.vcf.gz
    """

  stub:
    """
    touch ${sampleName}.vcf.gz    ${sampleName}.vcf.gz.tbi
    touch ${sampleName}.snps.vcf.gz   ${sampleName}.snps.vcf.gz.tbi
    touch ${sampleName}.indels.vcf.gz ${sampleName}.indels.vcf.gz.tbi
    """
}
```

- [ ] **Step 2: Update `makeSnpDensity` input tuple in `modules/cnv.nf`**

Find the `makeSnpDensity` process (line 242) and replace its `input:` block:

Old:
```nextflow
  input:
    tuple val(sampleName), path(freebayesVcfGz), path(freebayesVcfGzTbi), path(snpsVcfGz), path(snpsVcfGzTbi), path(indelsVcfGz), path(indelsVcfGzTbi), path(gvcfGz), path(gvcfGzTbi)
    tuple path(windows), path(genome)
```

New:
```nextflow
  input:
    tuple val(sampleName), path(vcfGz), path(vcfGzTbi), path(snpsVcfGz), path(snpsVcfGzTbi), path(indelsVcfGz), path(indelsVcfGzTbi)
    tuple path(windows), path(genome)
```

- [ ] **Step 3: Commit**

```bash
git add modules/snp.nf modules/cnv.nf
git commit -m "feat: add extractSampleVcf process; slim makeSnpDensity input tuple"
```

---

### Task 11: Rewrite workflow channel wiring in `processSingleExperiment.nf`

**Files:**
- Modify: `workflows/processSingleExperiment.nf`

- [ ] **Step 1: Update the `include` block at the top — remove dead processes, add new ones**

Replace the SNP include block (lines 20–27):

```nextflow
// SNP
include { makeRegionBed }                  from '../modules/snp.nf'
include { makeMultiSampleZeroCoverageBed } from '../modules/snp.nf'
include { freebayesMultiSample }           from '../modules/snp.nf'
include { splitGvcfAtZeroCoverage }        from '../modules/snp.nf'
include { concatMultiSampleVcf }           from '../modules/snp.nf'
include { makeConsensusFromGvcf }          from '../modules/snp.nf'
include { makeIndelTSV }                   from '../modules/snp.nf'
include { extractSampleVcf }              from '../modules/snp.nf'
```

- [ ] **Step 2: Replace the SNP section of the `ps` workflow body**

Remove everything from `freebayesResults = freebayes(...)` through the closing `freebayesMultiSample(...)` block (lines 90–128). Replace with:

```nextflow
    // ── Multi-sample FreeBayes (chunked) ─────────────────────────────────────

    regionsChannel = makeRegionBed(reorderFastaResults)
        .splitText()
        .map { it.trim() }
        .filter { it }

    allBams = gatkResults.bamTuple.map { sn, bam, bai -> bam }.collect()
    allBais = gatkResults.bamTuple.map { sn, bam, bai -> bai }.collect()

    // Per-region: joint multi-sample FreeBayes
    multiSampleChunks = freebayesMultiSample(allBams, allBais, reorderFastaResults, regionsChannel)

    // Per-region: union zero-coverage BED across all samples
    zeroCovBeds = makeMultiSampleZeroCoverageBed(allBams, allBais, regionsChannel)

    // Per-region: split gVCF at zero-coverage boundaries
    splitInput = multiSampleChunks.gvcf.join(zeroCovBeds, by: 0)
    splitChunks = splitGvcfAtZeroCoverage(splitInput, reorderFastaResults)

    // Gather all chunks → joint multi-sample VCF
    concatResult = concatMultiSampleVcf(
        splitChunks.map { regionKey, gvcf, tbi -> gvcf }.collect(),
        splitChunks.map { regionKey, gvcf, tbi -> tbi }.collect()
    )

    // Consensus FASTAs (one per sample)
    makeConsensusFromGvcf(concatResult.gvcf, reorderFastaResults)

    // Indels TSV (all samples from joint gVCF)
    makeIndelTSV(concatResult.gvcf)
        .collectFile(name: 'indels.tsv', storeDir: params.outputDir)

    // Per-sample VCF extraction for CNV steps
    sampleNames = gatkResults.bamTuple.map { sn, bam, bai -> sn }
    perSampleVcfResults = extractSampleVcf(
        sampleNames.combine(concatResult.vcf)
    )
```

- [ ] **Step 3: Update the CNV section to use `perSampleVcfResults` instead of `freebayesResults.vcf_files`**

Replace `makeSnpDensityResults = makeSnpDensity(freebayesResults.vcf_files, makeWindowFileResults)` with:

```nextflow
    makeSnpDensityResults = makeSnpDensity(perSampleVcfResults.vcf_files, makeWindowFileResults)
```

Replace the heterozygous SNP block (lines 157–170):

```nextflow
    if (params.ploidy != 1) {

      snpsVcf = perSampleVcfResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, snpsVcfGz, snpsVcfGzTbi)
      }

      convertedSnpsVcf = convertFreebayesToVarscanFormat(snpsVcf)

      getHeterozygousSNPsResults = getHeterozygousSNPs(convertedSnpsVcf)

      makeHeterozygousDensityBedResults = makeHeterozygousDensityBed(makeWindowFileResults, getHeterozygousSNPsResults)

      makeHeterozygousDensityBigwig(makeHeterozygousDensityBedResults, reorderFastaResults)
    }
```

- [ ] **Step 4: Verify parse with stub run**

```bash
nextflow run main.nf -entry processSingleExperiment -profile processSingleExperiment -stub 2>&1 | tail -30
```

Expected: processes listed, no `Caused by` parse errors. Stub workflow may error on channel mismatches — that's fine at this stage; watch for Groovy/DSL syntax errors specifically.

- [ ] **Step 5: Commit**

```bash
git add workflows/processSingleExperiment.nf
git commit -m "feat: rewire processSingleExperiment to use chunked multi-sample FreeBayes"
```

---

### Task 12: Remove dead processes from `modules/snp.nf`

**Files:**
- Modify: `modules/snp.nf`

- [ ] **Step 1: Delete the following process definitions from `modules/snp.nf`**

- `freebayes` (lines 4–66)
- `mergeVcfs` (lines 94–120)
- `makeMergedVariantIndex` (lines 122–148)
- `mergeGvcfs` (lines 182–214)
- `bcftoolsMpileupGvcf` (lines 150–180)

- [ ] **Step 2: Verify parse**

```bash
nextflow run main.nf -entry processSingleExperiment -profile processSingleExperiment -stub 2>&1 | grep -i "error\|No such process\|Missing" | head -20
```

Expected: no errors about missing processes (none of them are imported in the workflow after Task 11).

- [ ] **Step 3: Commit**

```bash
git add modules/snp.nf
git commit -m "chore: remove per-sample freebayes, mergeVcfs, makeMergedVariantIndex, mergeGvcfs, bcftoolsMpileupGvcf"
```

---

### Task 13: Full stub run validation

- [ ] **Step 1: Run the complete workflow in stub mode**

```bash
nextflow run main.nf -entry processSingleExperiment -profile processSingleExperiment -stub 2>&1 | tee /tmp/stub_run.log
```

- [ ] **Step 2: Verify all expected processes appear in the log**

```bash
grep "process\|Process\|\[" /tmp/stub_run.log | grep -v "^N E X T F L O W\|^Launching\|^executor\|^monitor" | head -40
```

Expected processes to appear: `makeRegionBed`, `makeMultiSampleZeroCoverageBed`, `freebayesMultiSample`, `splitGvcfAtZeroCoverage`, `concatMultiSampleVcf`, `makeConsensusFromGvcf`, `makeIndelTSV`, `extractSampleVcf`.

- [ ] **Step 3: Verify dead processes do NOT appear**

```bash
grep -E "freebayes\b|mergeVcfs|makeMergedVariantIndex|mergeGvcfs|bcftoolsMpileupGvcf" /tmp/stub_run.log
```

Expected: no output.

- [ ] **Step 4: Run the Python unit tests**

```bash
python -m pytest testing/t/test_splitGvcfAtZeroCoverage.py testing/t/test_makeConsensusFastaFromGvcf.py -v 2>&1 | tail -15
```

Expected: all tests PASSED.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "test: confirm full stub run and Python unit tests pass for multi-sample FreeBayes"
```
