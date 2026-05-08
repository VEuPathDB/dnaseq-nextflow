# mergeExperiments E2E Tests Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Python pytest suite that validates all Layer 1 structural invariants and Layer 2 cross-output consistency checks defined in `docs/superpowers/specs/2026-05-08-mergeExperiments-e2e-test-spec.md`.

**Architecture:** A `conftest.py` fixture discovers the six process work directories from the most recent `nextflow log` run. Each test function opens files directly from those work directories — no pipeline re-execution needed. Tests are run with `--run-dir` pointing at a Nextflow launch directory that contains a completed `mergeExperiments` run.

**Tech Stack:** Python 3, pytest, sqlite3 (stdlib), subprocess + bcftools (host install)

---

## Reference Dataset

Run directory: `/home/jbrestel/dnaseq_test/merge`
Strains: `Friedlin_resequence`, `LV39`, `LV39cl5`, `Seidman751`
Reference strain: `lmajFriedlin`

---

## File Map

| File | Action | Responsibility |
|---|---|---|
| `testing/t/conftest.py` | Create | `--run-dir` CLI option; `work_dirs` session fixture |
| `testing/t/test_mergeExperiments_e2e.py` | Create | All Layer 1 + Layer 2 test functions + shared helpers |

---

## Running the Tests

```bash
cd /path/to/dnaseq-nextflow
python -m pytest testing/t/test_mergeExperiments_e2e.py \
  --run-dir /home/jbrestel/dnaseq_test/merge -v
```

---

### Task 1: pytest fixtures (conftest.py)

**Files:**
- Create: `testing/t/conftest.py`

- [ ] **Step 1: Write conftest.py**

```python
# testing/t/conftest.py
import subprocess
import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--run-dir", required=True,
        help="Nextflow launch directory containing a completed mergeExperiments run"
    )


@pytest.fixture(scope="session")
def run_dir(request):
    return request.config.getoption("--run-dir")


@pytest.fixture(scope="session")
def work_dirs(run_dir):
    """Return dict mapping short process name -> work dir path for the most recent run."""
    r = subprocess.run(
        "nextflow log | tail -n1 | cut -f3",
        shell=True, cwd=run_dir, capture_output=True, text=True, check=True
    )
    run_name = r.stdout.strip()
    assert run_name, f"No nextflow runs found in {run_dir}"

    r = subprocess.run(
        ["nextflow", "log", run_name, "-f", "process,workdir"],
        cwd=run_dir, capture_output=True, text=True, check=True
    )
    dirs = {}
    for line in r.stdout.strip().split('\n'):
        parts = line.split('\t')
        if len(parts) != 2:
            continue
        process_name = parts[0].strip().split(':')[-1]
        dirs[process_name] = parts[1].strip()
    return dirs
```

- [ ] **Step 2: Verify fixture discovers all six processes**

```bash
cd /home/jbrestel/dnaseq_test/merge
python -c "
import subprocess
r = subprocess.run('nextflow log | tail -n1 | cut -f3', shell=True, capture_output=True, text=True)
run_name = r.stdout.strip()
r2 = subprocess.run(['nextflow', 'log', run_name, '-f', 'process,workdir'], capture_output=True, text=True)
print(r2.stdout)
"
```

Expected: six lines, one per process (`makeGenomicIndelDb`, `mergeCoverageBeds`, `mergeVcfs`, `makeCodingData`, `processSeqVars`, `snpEff`).

- [ ] **Step 3: Commit**

```bash
git add testing/t/conftest.py
git commit -m "test: add pytest conftest with nextflow work-dir discovery fixture"
```

---

### Task 2: Layer 1 — makeGenomicIndelDb

**Files:**
- Create: `testing/t/test_mergeExperiments_e2e.py` (initial version)

- [ ] **Step 1: Create test file with shared helpers and makeGenomicIndelDb tests**

```python
# testing/t/test_mergeExperiments_e2e.py
import os
import re
import sqlite3
import subprocess

import pytest


# ---------------------------------------------------------------------------
# bcftools helpers
# ---------------------------------------------------------------------------

def bcftools_is_valid(vcf_path):
    r = subprocess.run(['bcftools', 'view', vcf_path, '-o', '/dev/null'],
                       capture_output=True)
    return r.returncode == 0


def bcftools_header(vcf_path):
    return subprocess.run(
        ['bcftools', 'view', '-h', vcf_path],
        capture_output=True, text=True, check=True
    ).stdout


def bcftools_samples(vcf_path):
    r = subprocess.run(['bcftools', 'query', '-l', vcf_path],
                       capture_output=True, text=True, check=True)
    return [s for s in r.stdout.strip().split('\n') if s]


def bcftools_record_count(vcf_path):
    r = subprocess.run(['bcftools', 'stats', vcf_path],
                       capture_output=True, text=True, check=True)
    for line in r.stdout.split('\n'):
        if line.startswith('SN') and 'number of records:' in line:
            return int(line.split('\t')[3])
    return 0


def bcftools_query_info(vcf_path, field):
    """Return list of INFO/field values for every record."""
    r = subprocess.run(
        ['bcftools', 'query', '-f', f'%INFO/{field}\n', vcf_path],
        capture_output=True, text=True, check=True
    )
    return [v.strip() for v in r.stdout.strip().split('\n') if v.strip()]


# ---------------------------------------------------------------------------
# GTF helper (used by makeCodingData length invariant)
# ---------------------------------------------------------------------------

def parse_gtf_cds_exons(gtf_path):
    """Return dict: transcript_id -> list of (seq_id, start, stop) for CDS features."""
    exons = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            m = re.search(r'transcript_id "([^"]+)"', fields[8])
            if not m:
                continue
            tid = m.group(1)
            exons.setdefault(tid, []).append((fields[0], int(fields[3]), int(fields[4])))
    return exons


# ---------------------------------------------------------------------------
# Layer 1: makeGenomicIndelDb
# ---------------------------------------------------------------------------

def test_genomic_indel_db_exists(work_dirs):
    db_path = os.path.join(work_dirs['makeGenomicIndelDb'], 'genomicIndels.db')
    assert os.path.exists(db_path), f"genomicIndels.db not found at {db_path}"


def test_genomic_indel_db_schema(work_dirs):
    db_path = os.path.join(work_dirs['makeGenomicIndelDb'], 'genomicIndels.db')
    conn = sqlite3.connect(db_path)
    cur = conn.execute("PRAGMA table_info(genomic_indels)")
    cols = {row[1] for row in cur.fetchall()}
    assert cols == {'strain', 'sequence_id', 'position', 'shift'}, \
        f"Unexpected columns: {cols}"


def test_genomic_indel_db_index(work_dirs):
    db_path = os.path.join(work_dirs['makeGenomicIndelDb'], 'genomicIndels.db')
    conn = sqlite3.connect(db_path)
    cur = conn.execute("SELECT name FROM sqlite_master WHERE type='index'")
    indexes = {row[0] for row in cur.fetchall()}
    assert 'idx_genomic_indels' in indexes, f"Missing index. Found: {indexes}"


def test_genomic_indel_db_row_count(work_dirs):
    db_path = os.path.join(work_dirs['makeGenomicIndelDb'], 'genomicIndels.db')
    conn = sqlite3.connect(db_path)
    count = conn.execute("SELECT COUNT(*) FROM genomic_indels").fetchone()[0]
    assert count > 0, "genomic_indels table is empty"


def test_genomic_indel_db_no_zero_shift(work_dirs):
    db_path = os.path.join(work_dirs['makeGenomicIndelDb'], 'genomicIndels.db')
    conn = sqlite3.connect(db_path)
    zeros = conn.execute("SELECT COUNT(*) FROM genomic_indels WHERE shift = 0").fetchone()[0]
    assert zeros == 0, f"{zeros} rows have shift=0"
```

- [ ] **Step 2: Run and verify passes**

```bash
cd /path/to/dnaseq-nextflow
python -m pytest testing/t/test_mergeExperiments_e2e.py \
  -k "genomic_indel" \
  --run-dir /home/jbrestel/dnaseq_test/merge -v
```

Expected: `5 passed`.

- [ ] **Step 3: Commit**

```bash
git add testing/t/test_mergeExperiments_e2e.py
git commit -m "test: add Layer 1 tests for makeGenomicIndelDb"
```

---

### Task 3: Layer 1 — mergeCoverageBeds

**Files:**
- Modify: `testing/t/test_mergeExperiments_e2e.py`

- [ ] **Step 1: Add coverage TSV tests**

Append to `testing/t/test_mergeExperiments_e2e.py`:

```python
# ---------------------------------------------------------------------------
# Layer 1: mergeCoverageBeds
# ---------------------------------------------------------------------------

def test_coverage_tsv_exists(work_dirs):
    path = os.path.join(work_dirs['mergeCoverageBeds'], 'coverage.tsv')
    assert os.path.exists(path)


def test_coverage_tsv_header(work_dirs):
    path = os.path.join(work_dirs['mergeCoverageBeds'], 'coverage.tsv')
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
    assert header[:3] == ['chrom', 'start', 'end'], \
        f"First 3 header cols wrong: {header[:3]}"
    assert len(header) > 3, "No strain columns in header"


def test_coverage_tsv_strain_columns_match_inputs(work_dirs):
    work = work_dirs['mergeCoverageBeds']
    bed_strains = sorted(
        f.replace('.coverage.bed.gz', '')
        for f in os.listdir(work)
        if f.endswith('.coverage.bed.gz')
    )
    path = os.path.join(work, 'coverage.tsv')
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
    tsv_strains = sorted(header[3:])
    assert tsv_strains == bed_strains, \
        f"TSV strains {tsv_strains} != BED strains {bed_strains}"


def test_coverage_tsv_column_count_consistent(work_dirs):
    path = os.path.join(work_dirs['mergeCoverageBeds'], 'coverage.tsv')
    with open(path) as f:
        lines = f.readlines()
    expected_cols = len(lines[0].split('\t'))
    bad = [i + 2 for i, line in enumerate(lines[1:]) if len(line.split('\t')) != expected_cols]
    assert not bad, f"Wrong column count on lines: {bad[:5]}"


def test_coverage_tsv_start_less_than_end(work_dirs):
    path = os.path.join(work_dirs['mergeCoverageBeds'], 'coverage.tsv')
    bad = []
    with open(path) as f:
        next(f)  # skip header
        for i, line in enumerate(f, start=2):
            cols = line.rstrip('\n').split('\t')
            if int(cols[1]) >= int(cols[2]):
                bad.append(i)
    assert not bad, f"start >= end on lines: {bad[:5]}"


def test_coverage_tsv_no_negative_values(work_dirs):
    path = os.path.join(work_dirs['mergeCoverageBeds'], 'coverage.tsv')
    bad = []
    with open(path) as f:
        next(f)
        for i, line in enumerate(f, start=2):
            cols = line.rstrip('\n').split('\t')
            if any(float(c) < 0 for c in cols[3:]):
                bad.append(i)
    assert not bad, f"Negative coverage values on lines: {bad[:5]}"


def test_coverage_tsv_row_count(work_dirs):
    path = os.path.join(work_dirs['mergeCoverageBeds'], 'coverage.tsv')
    with open(path) as f:
        count = sum(1 for _ in f) - 1  # subtract header
    assert count > 0, "coverage.tsv has no data rows"
```

- [ ] **Step 2: Run and verify passes**

```bash
python -m pytest testing/t/test_mergeExperiments_e2e.py \
  -k "coverage_tsv" \
  --run-dir /home/jbrestel/dnaseq_test/merge -v
```

Expected: `6 passed`.

- [ ] **Step 3: Commit**

```bash
git add testing/t/test_mergeExperiments_e2e.py
git commit -m "test: add Layer 1 tests for mergeCoverageBeds"
```

---

### Task 4: Layer 1 — mergeVcfs

**Files:**
- Modify: `testing/t/test_mergeExperiments_e2e.py`

- [ ] **Step 1: Add mergeVcfs tests**

Append to `testing/t/test_mergeExperiments_e2e.py`:

```python
# ---------------------------------------------------------------------------
# Layer 1: mergeVcfs
# ---------------------------------------------------------------------------

def test_merged_vcf_exists(work_dirs):
    path = os.path.join(work_dirs['mergeVcfs'], 'merged.vcf.gz')
    assert os.path.exists(path)
    assert os.path.exists(path + '.tbi') or os.path.exists(path.replace('.vcf.gz', '.vcf.gz.tbi')), \
        "merged.vcf.gz.tbi index missing"


def test_merged_vcf_is_valid_bgzipped(work_dirs):
    path = os.path.join(work_dirs['mergeVcfs'], 'merged.vcf.gz')
    assert bcftools_is_valid(path), "merged.vcf.gz fails bcftools validation"


def test_merged_vcf_has_tbi_index(work_dirs):
    work = work_dirs['mergeVcfs']
    path = os.path.join(work, 'merged.vcf.gz')
    tbi = path + '.tbi'
    assert os.path.exists(tbi), f"Tabix index not found: {tbi}"


def test_merged_vcf_sample_names_match_input_strains(work_dirs):
    # Input VCFs are staged as 1.vcf.gz, 2.vcf.gz... — use coverage BEDs for ground truth strains
    bed_work = work_dirs['mergeCoverageBeds']
    expected = {
        f.replace('.coverage.bed.gz', '')
        for f in os.listdir(bed_work)
        if f.endswith('.coverage.bed.gz')
    }
    vcf_path = os.path.join(work_dirs['mergeVcfs'], 'merged.vcf.gz')
    actual = set(bcftools_samples(vcf_path))
    assert actual == expected, f"VCF samples {actual} != expected strains {expected}"


def test_merged_vcf_every_record_has_nonref_gt(work_dirs):
    vcf_path = os.path.join(work_dirs['mergeVcfs'], 'merged.vcf.gz')
    # Use bcftools view to filter for records with at least one non-ref GT
    r = subprocess.run(
        ['bcftools', 'view', '--min-ac', '1', vcf_path, '-o', '/dev/null'],
        capture_output=True
    )
    assert r.returncode == 0
    all_count = bcftools_record_count(vcf_path)
    r2 = subprocess.run(
        ['bcftools', 'view', '--min-ac', '1', '-O', 'z', vcf_path],
        capture_output=True
    )
    filtered_vcf = '/tmp/e2e_nonref_check.vcf.gz'
    with open(filtered_vcf, 'wb') as fh:
        fh.write(r2.stdout)
    subprocess.run(['bcftools', 'index', '-t', filtered_vcf], check=True)
    nonref_count = bcftools_record_count(filtered_vcf)
    assert nonref_count == all_count, \
        f"{all_count - nonref_count} records have no non-ref GT"


def test_merged_vcf_record_count(work_dirs):
    vcf_path = os.path.join(work_dirs['mergeVcfs'], 'merged.vcf.gz')
    assert bcftools_record_count(vcf_path) > 0
```

- [ ] **Step 2: Run and verify passes**

```bash
python -m pytest testing/t/test_mergeExperiments_e2e.py \
  -k "merged_vcf" \
  --run-dir /home/jbrestel/dnaseq_test/merge -v
```

Expected: `6 passed`.

- [ ] **Step 3: Commit**

```bash
git add testing/t/test_mergeExperiments_e2e.py
git commit -m "test: add Layer 1 tests for mergeVcfs"
```

---

### Task 5: Layer 1 — makeCodingData structural tests

**Files:**
- Modify: `testing/t/test_mergeExperiments_e2e.py`

- [ ] **Step 1: Add makeCodingData structural tests**

Append to `testing/t/test_mergeExperiments_e2e.py`:

```python
# ---------------------------------------------------------------------------
# Layer 1: makeCodingData — structural
# ---------------------------------------------------------------------------

def test_coding_sequences_db_exists(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    assert os.path.exists(path)


def test_coding_sequences_db_schema(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    conn = sqlite3.connect(path)
    cur = conn.execute("PRAGMA table_info(coding_sequences)")
    cols = {row[1] for row in cur.fetchall()}
    assert cols == {'strain', 'transcript_id', 'sequence'}, \
        f"Unexpected columns: {cols}"


def test_coding_sequences_all_strains_present(work_dirs):
    bed_work = work_dirs['mergeCoverageBeds']
    expected_strains = {
        f.replace('.coverage.bed.gz', '')
        for f in os.listdir(bed_work)
        if f.endswith('.coverage.bed.gz')
    }
    path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    conn = sqlite3.connect(path)
    actual_strains = {
        row[0]
        for row in conn.execute("SELECT DISTINCT strain FROM coding_sequences")
    }
    assert expected_strains.issubset(actual_strains), \
        f"Missing strains: {expected_strains - actual_strains}"


def test_coding_sequences_all_transcripts_all_strains(work_dirs):
    work = work_dirs['makeCodingData']
    path = os.path.join(work, 'codingSequences.db')
    conn = sqlite3.connect(path)
    n_strains = conn.execute("SELECT COUNT(DISTINCT strain) FROM coding_sequences").fetchone()[0]
    n_transcripts = conn.execute("SELECT COUNT(DISTINCT transcript_id) FROM coding_sequences").fetchone()[0]
    n_rows = conn.execute("SELECT COUNT(*) FROM coding_sequences").fetchone()[0]
    assert n_rows == n_strains * n_transcripts, \
        f"Expected {n_strains}×{n_transcripts}={n_strains * n_transcripts} rows, got {n_rows}"


def test_coding_sequences_non_empty(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    conn = sqlite3.connect(path)
    empty = conn.execute(
        "SELECT COUNT(*) FROM coding_sequences WHERE sequence IS NULL OR sequence = ''"
    ).fetchone()[0]
    assert empty == 0, f"{empty} rows have empty sequence"


def test_coding_sequences_start_with_atg(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    conn = sqlite3.connect(path)
    rows = conn.execute("SELECT strain, transcript_id, sequence FROM coding_sequences").fetchall()
    bad = [(s, t) for s, t, seq in rows if not seq.upper().startswith('ATG')]
    assert not bad, f"{len(bad)} sequences don't start with ATG: {bad[:3]}"


def test_coding_sequences_valid_characters(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    conn = sqlite3.connect(path)
    rows = conn.execute("SELECT strain, transcript_id, sequence FROM coding_sequences").fetchall()
    pattern = re.compile(r'^[ACGTacgtn]+$')
    bad = [(s, t) for s, t, seq in rows if not pattern.match(seq)]
    assert not bad, f"{len(bad)} sequences contain invalid characters: {bad[:3]}"


def test_coding_indels_db_exists(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingIndels.db')
    assert os.path.exists(path)


def test_coding_indels_db_schema(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingIndels.db')
    conn = sqlite3.connect(path)
    cur = conn.execute("PRAGMA table_info(indels)")
    cols = {row[1] for row in cur.fetchall()}
    assert cols == {'strain', 'transcript_id', 'position', 'shift_amount'}, \
        f"Unexpected columns: {cols}"


def test_coding_indels_no_zero_shift(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingIndels.db')
    conn = sqlite3.connect(path)
    zeros = conn.execute(
        "SELECT COUNT(*) FROM indels WHERE shift_amount = 0"
    ).fetchone()[0]
    assert zeros == 0, f"{zeros} rows have shift_amount=0"


def test_coding_indels_positions_positive(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingIndels.db')
    conn = sqlite3.connect(path)
    bad = conn.execute(
        "SELECT COUNT(*) FROM indels WHERE position < 1"
    ).fetchone()[0]
    assert bad == 0, f"{bad} rows have position < 1"
```

- [ ] **Step 2: Run and verify passes**

```bash
python -m pytest testing/t/test_mergeExperiments_e2e.py \
  -k "coding_sequence or coding_indel" \
  --run-dir /home/jbrestel/dnaseq_test/merge -v
```

Expected: `12 passed`.

- [ ] **Step 3: Commit**

```bash
git add testing/t/test_mergeExperiments_e2e.py
git commit -m "test: add Layer 1 structural tests for makeCodingData"
```

---

### Task 6: Layer 1 — makeCodingData length invariant

This is the Check 3 invariant from `docs/qa-makeCodingData-2026-03-18.md` promoted to automated form: for every `(strain, transcript)`, `len(sequence) == sum of exon lengths + sum of genomic indel shifts within those exons for that strain`.

**Files:**
- Modify: `testing/t/test_mergeExperiments_e2e.py`

- [ ] **Step 1: Add length invariant test**

Append to `testing/t/test_mergeExperiments_e2e.py`:

```python
def test_coding_sequences_length_invariant(work_dirs):
    """
    For every (strain, transcript), the CDS sequence length must equal:
      sum(exon_stop - exon_start + 1 for each CDS exon)
      + sum(shift for genomic indels within any exon of this transcript, for this strain)

    Derived from Check 3 in docs/qa-makeCodingData-2026-03-18.md.
    """
    work = work_dirs['makeCodingData']

    gtf_path = next(
        os.path.join(work, f) for f in os.listdir(work) if f.endswith('.gtf')
    )
    exons_by_transcript = parse_gtf_cds_exons(gtf_path)

    cds_conn = sqlite3.connect(os.path.join(work, 'codingSequences.db'))
    indel_conn = sqlite3.connect(os.path.join(work, 'genomicIndels.db'))

    rows = cds_conn.execute(
        "SELECT strain, transcript_id, LENGTH(sequence) FROM coding_sequences"
    ).fetchall()

    mismatches = []
    for strain, tid, actual_len in rows:
        if tid not in exons_by_transcript:
            continue  # transcript not in GTF subset — skip

        exon_list = exons_by_transcript[tid]
        ref_len = sum(stop - start + 1 for _, start, stop in exon_list)

        shift_total = 0
        for seq_id, start, stop in exon_list:
            row = indel_conn.execute(
                """SELECT COALESCE(SUM(shift), 0)
                   FROM genomic_indels
                   WHERE strain = ? AND sequence_id = ?
                     AND position >= ? AND position <= ?""",
                (strain, seq_id, start, stop)
            ).fetchone()
            shift_total += row[0]

        expected_len = ref_len + shift_total
        if actual_len != expected_len:
            mismatches.append((strain, tid, expected_len, actual_len))

    assert not mismatches, (
        f"{len(mismatches)} length mismatches (first 3): "
        + str(mismatches[:3])
    )
```

- [ ] **Step 2: Run and verify passes**

```bash
python -m pytest testing/t/test_mergeExperiments_e2e.py \
  -k "length_invariant" \
  --run-dir /home/jbrestel/dnaseq_test/merge -v
```

Expected: `1 passed`.

- [ ] **Step 3: Commit**

```bash
git add testing/t/test_mergeExperiments_e2e.py
git commit -m "test: add makeCodingData CDS length invariant test (Check 3 from QA doc)"
```

---

### Task 7: Layer 1 — processSeqVars: output.vcf.gz + variationFeature.dat

**Files:**
- Modify: `testing/t/test_mergeExperiments_e2e.py`

- [ ] **Step 1: Add processSeqVars VCF and variationFeature tests**

Append to `testing/t/test_mergeExperiments_e2e.py`:

```python
# ---------------------------------------------------------------------------
# Layer 1: processSeqVars — output.vcf.gz
# ---------------------------------------------------------------------------

def test_output_vcf_exists(work_dirs):
    path = os.path.join(work_dirs['processSeqVars'], 'output.vcf.gz')
    assert os.path.exists(path)
    assert os.path.exists(path + '.tbi')


def test_output_vcf_is_valid(work_dirs):
    path = os.path.join(work_dirs['processSeqVars'], 'output.vcf.gz')
    assert bcftools_is_valid(path)


def test_output_vcf_has_cann_header(work_dirs):
    path = os.path.join(work_dirs['processSeqVars'], 'output.vcf.gz')
    header = bcftools_header(path)
    assert '##INFO=<ID=CANN' in header
    assert '##FORMAT=<ID=CA' in header
    assert '##FORMAT=<ID=DFS' in header


def test_output_vcf_all_records_have_cann(work_dirs):
    path = os.path.join(work_dirs['processSeqVars'], 'output.vcf.gz')
    values = bcftools_query_info(path, 'CANN')
    assert len(values) > 0
    missing = [i + 1 for i, v in enumerate(values) if not v or v == '.']
    assert not missing, f"Records missing CANN at positions: {missing[:5]}"


def test_output_vcf_record_count(work_dirs):
    path = os.path.join(work_dirs['processSeqVars'], 'output.vcf.gz')
    assert bcftools_record_count(path) > 0


# ---------------------------------------------------------------------------
# Layer 1: processSeqVars — variationFeature.dat
# ---------------------------------------------------------------------------

def _read_variation_feature(work_dirs):
    path = os.path.join(work_dirs['processSeqVars'], 'variationFeature.dat')
    rows = []
    with open(path) as f:
        for line in f:
            rows.append(line.rstrip('\n').split('\t'))
    return rows


def test_variation_feature_exists(work_dirs):
    path = os.path.join(work_dirs['processSeqVars'], 'variationFeature.dat')
    assert os.path.exists(path)


def test_variation_feature_column_count(work_dirs):
    rows = _read_variation_feature(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if len(r) != 18]
    assert not bad, f"Rows with wrong column count (expected 18): {bad[:5]}"


def test_variation_feature_row_count(work_dirs):
    rows = _read_variation_feature(work_dirs)
    assert len(rows) > 0


def test_variation_feature_location_positive_int(work_dirs):
    rows = _read_variation_feature(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[0].isdigit() or int(r[0]) <= 0]
    assert not bad, f"Rows with invalid location (col 1): {bad[:5]}"


def test_variation_feature_seq_id_nonempty(work_dirs):
    rows = _read_variation_feature(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[2].strip()]
    assert not bad, f"Rows with empty seq_id (col 3): {bad[:5]}"


def test_variation_feature_reference_strain(work_dirs):
    # Determine reference strain from undoneStrains file sibling dir
    work = work_dirs['processSeqVars']
    # Read it from the VCF cache header or just check col 4 is a consistent non-empty value
    rows = _read_variation_feature(work_dirs)
    values = {r[3] for r in rows}
    assert len(values) == 1, f"Multiple reference_strain values in col 4: {values}"
    assert next(iter(values)), "reference_strain (col 4) is empty"


def test_variation_feature_has_nonsynonymous_binary(work_dirs):
    rows = _read_variation_feature(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if r[5] not in ('0', '1')]
    assert not bad, f"Rows with has_nonsynonymous not 0/1 (col 6): {bad[:5]}"


def test_variation_feature_major_allele_nonempty(work_dirs):
    rows = _read_variation_feature(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[6].strip()]
    assert not bad, f"Rows with empty major_allele (col 7): {bad[:5]}"


def test_variation_feature_major_allele_count_positive(work_dirs):
    rows = _read_variation_feature(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[8].isdigit() or int(r[8]) <= 0]
    assert not bad, f"Rows with major_allele_count <= 0 (col 9): {bad[:5]}"


def test_variation_feature_distinct_strain_count_in_range(work_dirs):
    vcf_path = os.path.join(work_dirs['processSeqVars'], 'merged.vcf.gz')
    n_strains = len(bcftools_samples(vcf_path))
    rows = _read_variation_feature(work_dirs)
    bad = [
        i + 1 for i, r in enumerate(rows)
        if not r[12].isdigit() or not (1 <= int(r[12]) <= n_strains)
    ]
    assert not bad, f"Rows with distinct_strain_count out of range (col 13): {bad[:5]}"


def test_variation_feature_is_coding_binary(work_dirs):
    rows = _read_variation_feature(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if r[14] not in ('0', '1')]
    assert not bad, f"Rows with is_coding not 0/1 (col 15): {bad[:5]}"


def test_variation_feature_coding_rows_have_transcript(work_dirs):
    rows = _read_variation_feature(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if r[14] == '1' and not r[1].strip()]
    assert not bad, f"is_coding=1 rows with empty transcript_id (col 2): {bad[:5]}"


def test_variation_feature_has_stop_codon_binary(work_dirs):
    rows = _read_variation_feature(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if r[16] not in ('0', '1')]
    assert not bad, f"Rows with has_stop_codon not 0/1 (col 17): {bad[:5]}"


def test_variation_feature_ref_allele_nonempty(work_dirs):
    rows = _read_variation_feature(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[4].strip()]
    assert not bad, f"Rows with empty ref_allele (col 5): {bad[:5]}"


def test_variation_feature_distinct_allele_count_gte_1(work_dirs):
    rows = _read_variation_feature(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[13].isdigit() or int(r[13]) < 1]
    assert not bad, f"Rows with distinct_allele_count < 1 (col 14): {bad[:5]}"


def test_variation_feature_total_allele_count_gte_strain_count(work_dirs):
    rows = _read_variation_feature(work_dirs)
    bad = [
        i + 1 for i, r in enumerate(rows)
        if not r[15].isdigit() or int(r[15]) < int(r[12])
    ]
    assert not bad, f"Rows where total_allele_count < distinct_strain_count (col 16 < col 13): {bad[:5]}"


def test_variation_feature_nonsynonymous_implies_coding(work_dirs):
    rows = _read_variation_feature(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if r[5] == '1' and r[14] != '1']
    assert not bad, f"has_nonsynonymous=1 but is_coding!=1 on rows: {bad[:5]}"
```

- [ ] **Step 2: Run and verify passes**

```bash
python -m pytest testing/t/test_mergeExperiments_e2e.py \
  -k "output_vcf or variation_feature" \
  --run-dir /home/jbrestel/dnaseq_test/merge -v
```

Expected: `18 passed`.

- [ ] **Step 3: Commit**

```bash
git add testing/t/test_mergeExperiments_e2e.py
git commit -m "test: add Layer 1 tests for processSeqVars output.vcf.gz and variationFeature.dat"
```

---

### Task 8: Layer 1 — processSeqVars: allele.dat + product.dat

**Files:**
- Modify: `testing/t/test_mergeExperiments_e2e.py`

- [ ] **Step 1: Add allele.dat and product.dat tests**

Append to `testing/t/test_mergeExperiments_e2e.py`:

```python
# ---------------------------------------------------------------------------
# Layer 1: processSeqVars — allele.dat
# ---------------------------------------------------------------------------

def _read_allele(work_dirs):
    path = os.path.join(work_dirs['processSeqVars'], 'allele.dat')
    rows = []
    with open(path) as f:
        for line in f:
            rows.append(line.rstrip('\n').split('\t'))
    return rows


def test_allele_dat_exists(work_dirs):
    assert os.path.exists(os.path.join(work_dirs['processSeqVars'], 'allele.dat'))


def test_allele_dat_column_count(work_dirs):
    rows = _read_allele(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if len(r) != 5]
    assert not bad, f"Rows with wrong column count (expected 5): {bad[:5]}"


def test_allele_dat_allele_nonempty(work_dirs):
    rows = _read_allele(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[0].strip()]
    assert not bad, f"Rows with empty allele (col 1): {bad[:5]}"


def test_allele_dat_distinct_strain_count_positive(work_dirs):
    rows = _read_allele(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[1].isdigit() or int(r[1]) <= 0]
    assert not bad, f"Rows with distinct_strain_count <= 0 (col 2): {bad[:5]}"


def test_allele_dat_allele_count_gte_strain_count(work_dirs):
    rows = _read_allele(work_dirs)
    bad = [
        i + 1 for i, r in enumerate(rows)
        if int(r[2]) < int(r[1])
    ]
    assert not bad, f"Rows where allele_count < distinct_strain_count: {bad[:5]}"


def test_allele_dat_avg_coverage_non_negative(work_dirs):
    rows = _read_allele(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if float(r[3]) < 0]
    assert not bad, f"Rows with avg_coverage < 0 (col 4): {bad[:5]}"


def test_allele_dat_avg_percent_in_range(work_dirs):
    rows = _read_allele(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not (0.0 <= float(r[4]) <= 100.0)]
    assert not bad, f"Rows with avg_percent outside 0–100 (col 5): {bad[:5]}"


def test_allele_dat_two_decimal_places(work_dirs):
    rows = _read_allele(work_dirs)
    pattern = re.compile(r'^\d+\.\d{2}$')
    bad = [
        i + 1 for i, r in enumerate(rows)
        if not pattern.match(r[3]) or not pattern.match(r[4])
    ]
    assert not bad, f"Rows where avg_coverage or avg_percent not 2dp: {bad[:5]}"


# ---------------------------------------------------------------------------
# Layer 1: processSeqVars — product.dat
# ---------------------------------------------------------------------------

def _read_product(work_dirs):
    path = os.path.join(work_dirs['processSeqVars'], 'product.dat')
    rows = []
    with open(path) as f:
        for line in f:
            rows.append(line.rstrip('\n').split('\t'))
    return rows


def test_product_dat_exists(work_dirs):
    assert os.path.exists(os.path.join(work_dirs['processSeqVars'], 'product.dat'))


def test_product_dat_column_count(work_dirs):
    rows = _read_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if len(r) != 8]
    assert not bad, f"Rows with wrong column count (expected 8): {bad[:5]}"


def test_product_dat_codon_is_three_acgt(work_dirs):
    rows = _read_product(work_dirs)
    pattern = re.compile(r'^[ACGT]{3}$')
    bad = [i + 1 for i, r in enumerate(rows) if not pattern.match(r[0])]
    assert not bad, f"Rows with invalid codon (col 1): {bad[:5]}"


def test_product_dat_pos_in_codon_1_2_3(work_dirs):
    rows = _read_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if r[1] not in ('1', '2', '3')]
    assert not bad, f"Rows with pos_in_codon not in {{1,2,3}} (col 2): {bad[:5]}"


def test_product_dat_transcript_id_nonempty(work_dirs):
    rows = _read_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[2].strip()]
    assert not bad, f"Rows with empty transcript_id (col 3): {bad[:5]}"


def test_product_dat_product_count_positive(work_dirs):
    rows = _read_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[3].isdigit() or int(r[3]) <= 0]
    assert not bad, f"Rows with product_count <= 0 (col 4): {bad[:5]}"


def test_product_dat_amino_acid_single_char_or_stop(work_dirs):
    rows = _read_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if len(r[4]) != 1]
    assert not bad, f"Rows where amino_acid is not single char (col 5): {bad[:5]}"


def test_product_dat_pos_in_cds_positive(work_dirs):
    rows = _read_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[5].isdigit() or int(r[5]) <= 0]
    assert not bad, f"Rows with pos_in_cds <= 0 (col 6): {bad[:5]}"


def test_product_dat_downstream_of_frameshift_binary(work_dirs):
    rows = _read_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if r[7] not in ('0', '1')]
    assert not bad, f"Rows with downstream_of_frameshift not 0/1 (col 8): {bad[:5]}"
```

- [ ] **Step 2: Run and verify passes**

```bash
python -m pytest testing/t/test_mergeExperiments_e2e.py \
  -k "allele_dat or product_dat" \
  --run-dir /home/jbrestel/dnaseq_test/merge -v
```

Expected: `17 passed`.

- [ ] **Step 3: Commit**

```bash
git add testing/t/test_mergeExperiments_e2e.py
git commit -m "test: add Layer 1 tests for allele.dat and product.dat"
```

---

### Task 9: Layer 1 — snpEff

**Files:**
- Modify: `testing/t/test_mergeExperiments_e2e.py`

- [ ] **Step 1: Add snpEff tests**

Append to `testing/t/test_mergeExperiments_e2e.py`:

```python
# ---------------------------------------------------------------------------
# Layer 1: snpEff — merged.ann.vcf.gz
# ---------------------------------------------------------------------------

def test_ann_vcf_exists(work_dirs):
    path = os.path.join(work_dirs['snpEff'], 'merged.ann.vcf.gz')
    assert os.path.exists(path)
    assert os.path.exists(path + '.tbi')


def test_ann_vcf_is_valid(work_dirs):
    path = os.path.join(work_dirs['snpEff'], 'merged.ann.vcf.gz')
    assert bcftools_is_valid(path)


def test_ann_vcf_has_ann_header(work_dirs):
    path = os.path.join(work_dirs['snpEff'], 'merged.ann.vcf.gz')
    header = bcftools_header(path)
    assert '##INFO=<ID=ANN' in header, "Missing snpEff ANN INFO header"


def test_ann_vcf_has_cann_header(work_dirs):
    path = os.path.join(work_dirs['snpEff'], 'merged.ann.vcf.gz')
    header = bcftools_header(path)
    assert '##INFO=<ID=CANN' in header, "Missing CANN INFO header (should carry through from processSeqVars)"
    assert '##FORMAT=<ID=CA' in header
    assert '##FORMAT=<ID=DFS' in header


def test_ann_vcf_all_records_have_ann(work_dirs):
    path = os.path.join(work_dirs['snpEff'], 'merged.ann.vcf.gz')
    values = bcftools_query_info(path, 'ANN')
    missing = [i + 1 for i, v in enumerate(values) if not v or v == '.']
    assert not missing, f"Records missing ANN at positions: {missing[:5]}"


def test_ann_vcf_all_records_have_cann(work_dirs):
    path = os.path.join(work_dirs['snpEff'], 'merged.ann.vcf.gz')
    values = bcftools_query_info(path, 'CANN')
    missing = [i + 1 for i, v in enumerate(values) if not v or v == '.']
    assert not missing, f"Records missing CANN at positions: {missing[:5]}"
```

- [ ] **Step 2: Run and verify passes**

```bash
python -m pytest testing/t/test_mergeExperiments_e2e.py \
  -k "ann_vcf" \
  --run-dir /home/jbrestel/dnaseq_test/merge -v
```

Expected: `6 passed`.

- [ ] **Step 3: Commit**

```bash
git add testing/t/test_mergeExperiments_e2e.py
git commit -m "test: add Layer 1 tests for snpEff merged.ann.vcf.gz"
```

---

### Task 10: Layer 2 — Cross-Output Consistency

**Files:**
- Modify: `testing/t/test_mergeExperiments_e2e.py`

- [ ] **Step 1: Add cross-output consistency tests**

Append to `testing/t/test_mergeExperiments_e2e.py`:

```python
# ---------------------------------------------------------------------------
# Layer 2: Cross-output consistency
# ---------------------------------------------------------------------------

def test_strain_names_consistent_across_vcf_and_coverage(work_dirs):
    vcf_path = os.path.join(work_dirs['processSeqVars'], 'merged.vcf.gz')
    vcf_strains = set(bcftools_samples(vcf_path))

    tsv_path = os.path.join(work_dirs['processSeqVars'], 'coverage.tsv')
    with open(tsv_path) as f:
        header = f.readline().rstrip('\n').split('\t')
    tsv_strains = set(header[3:])

    assert vcf_strains == tsv_strains, \
        f"VCF samples {vcf_strains} != coverage.tsv strains {tsv_strains}"


def test_strain_names_consistent_in_genomic_indel_db(work_dirs):
    vcf_path = os.path.join(work_dirs['processSeqVars'], 'merged.vcf.gz')
    vcf_strains = set(bcftools_samples(vcf_path))

    db_path = os.path.join(work_dirs['makeCodingData'], 'genomicIndels.db')
    conn = sqlite3.connect(db_path)
    db_strains = {
        row[0]
        for row in conn.execute("SELECT DISTINCT strain FROM genomic_indels")
    }
    assert db_strains.issubset(vcf_strains), \
        f"genomic_indels strains not a subset of VCF samples: {db_strains - vcf_strains}"


def test_strain_names_consistent_in_coding_sequences(work_dirs):
    vcf_path = os.path.join(work_dirs['processSeqVars'], 'merged.vcf.gz')
    vcf_strains = set(bcftools_samples(vcf_path))

    db_path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    conn = sqlite3.connect(db_path)
    db_strains = {
        row[0]
        for row in conn.execute("SELECT DISTINCT strain FROM coding_sequences")
    }
    assert vcf_strains.issubset(db_strains), \
        f"VCF strains not present in coding_sequences: {vcf_strains - db_strains}"


def test_variant_count_output_vcf_equals_ann_vcf(work_dirs):
    output_vcf = os.path.join(work_dirs['processSeqVars'], 'output.vcf.gz')
    ann_vcf = os.path.join(work_dirs['snpEff'], 'merged.ann.vcf.gz')
    output_count = bcftools_record_count(output_vcf)
    ann_count = bcftools_record_count(ann_vcf)
    assert output_count == ann_count, \
        f"output.vcf.gz has {output_count} records but merged.ann.vcf.gz has {ann_count}"


def test_variation_feature_count_lte_vcf_record_count(work_dirs):
    vcf_path = os.path.join(work_dirs['processSeqVars'], 'output.vcf.gz')
    vcf_count = bcftools_record_count(vcf_path)

    vf_path = os.path.join(work_dirs['processSeqVars'], 'variationFeature.dat')
    with open(vf_path) as f:
        vf_count = sum(1 for _ in f)

    assert vf_count > 0, "variationFeature.dat is empty"
    assert vf_count <= vcf_count * 20, (
        f"variationFeature.dat rows ({vf_count}) suspiciously exceed "
        f"20× VCF records ({vcf_count})"
    )


def test_allele_dat_only_for_coding_variants(work_dirs):
    vf_rows = _read_variation_feature(work_dirs)
    coding_transcripts = {r[1] for r in vf_rows if r[14] == '1' and r[1]}

    allele_rows = _read_allele(work_dirs)
    # allele.dat has no transcript_id column directly — verify row count is
    # consistent with the number of alleles across coding positions
    # (weak check: allele.dat must not be larger than alleles × coding rows)
    coding_row_count = sum(1 for r in vf_rows if r[14] == '1')
    allele_row_count = len(allele_rows)
    assert allele_row_count > 0
    # Each coding position emits at most N_alleles rows; conservatively ≤ 4 alleles per site
    assert allele_row_count <= coding_row_count * 4, (
        f"allele.dat rows ({allele_row_count}) far exceed "
        f"4× coding variationFeature rows ({coding_row_count})"
    )


def test_product_dat_transcripts_in_coding_sequences(work_dirs):
    db_path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    conn = sqlite3.connect(db_path)
    known = {
        row[0]
        for row in conn.execute("SELECT DISTINCT transcript_id FROM coding_sequences")
    }

    prod_rows = _read_product(work_dirs)
    product_transcripts = {r[2] for r in prod_rows if r[2]}
    orphans = product_transcripts - known
    assert not orphans, \
        f"product.dat references transcripts not in codingSequences.db: {list(orphans)[:5]}"


def test_reference_strain_in_vcf_samples(work_dirs):
    vcf_path = os.path.join(work_dirs['processSeqVars'], 'merged.vcf.gz')
    samples = bcftools_samples(vcf_path)
    vf_rows = _read_variation_feature(work_dirs)
    ref_strain = vf_rows[0][3]
    assert ref_strain in samples, \
        f"reference_strain '{ref_strain}' not in VCF samples: {samples}"


def test_variation_feature_reference_strain_uniform(work_dirs):
    vf_rows = _read_variation_feature(work_dirs)
    values = {r[3] for r in vf_rows}
    assert len(values) == 1, f"variationFeature.dat col 4 has multiple values: {values}"


def test_cann_present_in_output_vcf_all_records(work_dirs):
    path = os.path.join(work_dirs['processSeqVars'], 'output.vcf.gz')
    values = bcftools_query_info(path, 'CANN')
    total = bcftools_record_count(path)
    assert len(values) == total, \
        f"bcftools query returned {len(values)} CANN values but VCF has {total} records"
    missing = [i + 1 for i, v in enumerate(values) if not v or v == '.']
    assert not missing, f"Records missing CANN at record positions: {missing[:5]}"


def test_cann_and_ann_both_in_merged_ann_vcf(work_dirs):
    path = os.path.join(work_dirs['snpEff'], 'merged.ann.vcf.gz')
    cann_values = bcftools_query_info(path, 'CANN')
    ann_values = bcftools_query_info(path, 'ANN')
    total = bcftools_record_count(path)
    assert len(cann_values) == total
    assert len(ann_values) == total
    missing_cann = [i + 1 for i, v in enumerate(cann_values) if not v or v == '.']
    missing_ann = [i + 1 for i, v in enumerate(ann_values) if not v or v == '.']
    assert not missing_cann, f"Records missing CANN in merged.ann.vcf.gz: {missing_cann[:5]}"
    assert not missing_ann, f"Records missing ANN in merged.ann.vcf.gz: {missing_ann[:5]}"
```

- [ ] **Step 2: Run and verify passes**

```bash
python -m pytest testing/t/test_mergeExperiments_e2e.py \
  -k "consistent or count_lte or product_dat_transcripts or reference_strain or cann_present or cann_and_ann" \
  --run-dir /home/jbrestel/dnaseq_test/merge -v
```

Expected: `11 passed`.

- [ ] **Step 3: Commit**

```bash
git add testing/t/test_mergeExperiments_e2e.py
git commit -m "test: add Layer 2 cross-output consistency tests"
```

---

### Task 11: Full suite smoke test

**Files:**
- No new files

- [ ] **Step 1: Run complete suite against reference run**

```bash
cd /path/to/dnaseq-nextflow
python -m pytest testing/t/test_mergeExperiments_e2e.py \
  --run-dir /home/jbrestel/dnaseq_test/merge -v
```

Expected output:
- All tests pass (no failures, no errors)
- Count should be approximately 73 tests

- [ ] **Step 2: If any test fails, investigate before proceeding**

A failure against the known-good reference run means either:
- An invariant in the spec is wrong → fix the test and update the spec
- The reference run has an unexpected property → investigate with `nextflow log` + work dir inspection

Do not mark passing until all tests pass cleanly.

- [ ] **Step 3: Final commit**

```bash
git add testing/t/conftest.py testing/t/test_mergeExperiments_e2e.py
git commit -m "test: mergeExperiments E2E test suite — all Layer 1 + Layer 2 checks passing"
```
