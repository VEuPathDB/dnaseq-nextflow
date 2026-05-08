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
