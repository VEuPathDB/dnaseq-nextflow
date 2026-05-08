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
    raise ValueError(f"Could not parse record count from bcftools stats output for {vcf_path}")


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
    with sqlite3.connect(db_path) as conn:
        cur = conn.execute("PRAGMA table_info(genomic_indels)")
        cols = {row[1] for row in cur.fetchall()}
    assert cols == {'strain', 'sequence_id', 'position', 'shift'}, \
        f"Unexpected columns: {cols}"


def test_genomic_indel_db_index(work_dirs):
    db_path = os.path.join(work_dirs['makeGenomicIndelDb'], 'genomicIndels.db')
    with sqlite3.connect(db_path) as conn:
        cur = conn.execute("SELECT name FROM sqlite_master WHERE type='index'")
        indexes = {row[0] for row in cur.fetchall()}
    assert 'idx_genomic_indels' in indexes, f"Missing index. Found: {indexes}"


def test_genomic_indel_db_row_count(work_dirs):
    db_path = os.path.join(work_dirs['makeGenomicIndelDb'], 'genomicIndels.db')
    with sqlite3.connect(db_path) as conn:
        count = conn.execute("SELECT COUNT(*) FROM genomic_indels").fetchone()[0]
    assert count > 0, "genomic_indels table is empty"


def test_genomic_indel_db_no_zero_shift(work_dirs):
    db_path = os.path.join(work_dirs['makeGenomicIndelDb'], 'genomicIndels.db')
    with sqlite3.connect(db_path) as conn:
        zeros = conn.execute("SELECT COUNT(*) FROM genomic_indels WHERE shift = 0").fetchone()[0]
    assert zeros == 0, f"{zeros} rows have shift=0"


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
