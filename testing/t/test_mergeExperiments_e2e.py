import os
import re
import sqlite3
import subprocess
import tempfile

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


# ---------------------------------------------------------------------------
# Layer 1: mergeVcfs
# ---------------------------------------------------------------------------

def test_merged_vcf_exists(work_dirs):
    path = os.path.join(work_dirs['mergeVcfs'], 'merged.vcf.gz')
    assert os.path.exists(path)


def test_merged_vcf_is_valid_bgzipped(work_dirs):
    path = os.path.join(work_dirs['mergeVcfs'], 'merged.vcf.gz')
    assert bcftools_is_valid(path), "merged.vcf.gz fails bcftools validation"


def test_merged_vcf_input_vcfs_have_tbi_index(work_dirs):
    """The merge process indexes each input VCF before merging."""
    work = work_dirs['mergeVcfs']
    input_vcfs = sorted(f for f in os.listdir(work) if f != 'merged.vcf.gz' and f.endswith('.vcf.gz'))
    missing = [v for v in input_vcfs if not os.path.exists(os.path.join(work, v + '.tbi'))]
    assert not missing, f"Input VCFs missing tabix index: {missing}"


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


def test_merged_vcf_mostly_nonref_gt(work_dirs):
    """Nearly all records should have at least one non-ref GT.
    Allow ≤1% for records where all samples are ./. or 0/0 (valid FreeBayes output
    at low-coverage sites)."""
    vcf_path = os.path.join(work_dirs['mergeVcfs'], 'merged.vcf.gz')
    all_count = bcftools_record_count(vcf_path)
    with tempfile.NamedTemporaryFile(suffix='.vcf.gz', delete=False) as tmp:
        filtered_vcf = tmp.name
    try:
        subprocess.run(
            ['bcftools', 'view', '--min-ac', '1', '-O', 'z', '-o', filtered_vcf, vcf_path],
            capture_output=True, check=True
        )
        subprocess.run(['bcftools', 'index', '-t', filtered_vcf], check=True)
        nonref_count = bcftools_record_count(filtered_vcf)
    finally:
        for f in [filtered_vcf, filtered_vcf + '.tbi']:
            if os.path.exists(f):
                os.unlink(f)
    tolerance = int(all_count * 0.01)
    assert all_count - nonref_count <= tolerance, \
        f"{all_count - nonref_count} records have no non-ref GT (>{tolerance} allowed, {all_count} total)"


def test_merged_vcf_record_count(work_dirs):
    vcf_path = os.path.join(work_dirs['mergeVcfs'], 'merged.vcf.gz')
    assert bcftools_record_count(vcf_path) > 0


# ---------------------------------------------------------------------------
# Layer 1: makeCodingData — structural
# ---------------------------------------------------------------------------

def test_coding_sequences_db_exists(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    assert os.path.exists(path)


def test_coding_sequences_db_schema(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    with sqlite3.connect(path) as conn:
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
    with sqlite3.connect(path) as conn:
        actual_strains = {
            row[0]
            for row in conn.execute("SELECT DISTINCT strain FROM coding_sequences")
        }
    assert expected_strains.issubset(actual_strains), \
        f"Missing strains: {expected_strains - actual_strains}"


def test_coding_sequences_all_transcripts_all_strains(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    with sqlite3.connect(path) as conn:
        n_strains = conn.execute("SELECT COUNT(DISTINCT strain) FROM coding_sequences").fetchone()[0]
        n_transcripts = conn.execute("SELECT COUNT(DISTINCT transcript_id) FROM coding_sequences").fetchone()[0]
        n_rows = conn.execute("SELECT COUNT(*) FROM coding_sequences").fetchone()[0]
    assert n_rows == n_strains * n_transcripts, \
        f"Expected {n_strains}×{n_transcripts}={n_strains * n_transcripts} rows, got {n_rows}"


def test_coding_sequences_non_empty(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    with sqlite3.connect(path) as conn:
        empty = conn.execute(
            "SELECT COUNT(*) FROM coding_sequences WHERE sequence IS NULL OR sequence = ''"
        ).fetchone()[0]
    assert empty == 0, f"{empty} rows have empty sequence"


def test_coding_sequences_start_with_atg(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    with sqlite3.connect(path) as conn:
        rows = conn.execute("SELECT strain, transcript_id, sequence FROM coding_sequences").fetchall()
    bad = [(s, t) for s, t, seq in rows if not seq.upper().startswith('ATG')]
    assert not bad, f"{len(bad)} sequences don't start with ATG: {bad[:3]}"


def test_coding_sequences_valid_characters(work_dirs):
    """Sequences may contain standard IUPAC nucleotide characters including:
    - ACGT: standard bases
    - N: masked/unknown base
    - X: het indel placeholder (non-IUPAC, pipeline-specific)
    - RYSWKM: IUPAC ambiguity codes for heterozygous SNP positions
    """
    path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    with sqlite3.connect(path) as conn:
        rows = conn.execute("SELECT strain, transcript_id, sequence FROM coding_sequences").fetchall()
    pattern = re.compile(r'^[ACGTRYSWKMNXacgtryswkmnx]+$')
    bad = [(s, t) for s, t, seq in rows if not pattern.match(seq)]
    assert not bad, f"{len(bad)} sequences contain invalid characters: {bad[:3]}"


def test_coding_indels_db_exists(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingIndels.db')
    assert os.path.exists(path)


def test_coding_indels_db_schema(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingIndels.db')
    with sqlite3.connect(path) as conn:
        cur = conn.execute("PRAGMA table_info(indels)")
        cols = {row[1] for row in cur.fetchall()}
    assert cols == {'strain', 'transcript_id', 'position', 'shift_amount'}, \
        f"Unexpected columns: {cols}"


def test_coding_indels_no_zero_shift(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingIndels.db')
    with sqlite3.connect(path) as conn:
        zeros = conn.execute(
            "SELECT COUNT(*) FROM indels WHERE shift_amount = 0"
        ).fetchone()[0]
    assert zeros == 0, f"{zeros} rows have shift_amount=0"


def test_coding_indels_positions_positive(work_dirs):
    path = os.path.join(work_dirs['makeCodingData'], 'codingIndels.db')
    with sqlite3.connect(path) as conn:
        bad = conn.execute(
            "SELECT COUNT(*) FROM indels WHERE position < 1"
        ).fetchone()[0]
    assert bad == 0, f"{bad} rows have position < 1"


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

    with sqlite3.connect(os.path.join(work, 'codingSequences.db')) as cds_conn, \
         sqlite3.connect(os.path.join(work, 'genomicIndels.db')) as indel_conn:
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


def test_output_vcf_has_cann_header_and_some_values(work_dirs):
    """CANN is only set on records overlapping coding regions.
    Non-coding variants have CANN='.'. Verify: (1) CANN header present,
    (2) at least some records have non-'.' CANN (i.e., the pipeline did produce
    coding annotations)."""
    path = os.path.join(work_dirs['processSeqVars'], 'output.vcf.gz')
    header = bcftools_header(path)
    assert '##INFO=<ID=CANN' in header, "CANN INFO header missing"
    values = bcftools_query_info(path, 'CANN')
    annotated = [v for v in values if v and v != '.']
    assert len(annotated) > 0, "No records have CANN annotation (expected at least some coding variants)"


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
    """distinct_strain_count (col 13) includes the reference strain as a separate variation,
    so the upper bound is N_vcf_samples + 1."""
    vcf_path = os.path.join(work_dirs['processSeqVars'], 'merged.vcf.gz')
    n_strains = len(bcftools_samples(vcf_path))
    rows = _read_variation_feature(work_dirs)
    bad = [
        i + 1 for i, r in enumerate(rows)
        if not r[12].isdigit() or not (1 <= int(r[12]) <= n_strains + 1)
    ]
    assert not bad, f"Rows with distinct_strain_count out of range [1,{n_strains + 1}] (col 13): {bad[:5]}"


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


def test_allele_dat_row_count(work_dirs):
    rows = _read_allele(work_dirs)
    assert len(rows) > 0, "allele.dat has no data rows"


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
        if not r[1].isdigit() or not r[2].isdigit() or int(r[2]) < int(r[1])
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


def test_product_dat_row_count(work_dirs):
    rows = _read_product(work_dirs)
    assert len(rows) > 0, "product.dat has no data rows"


def test_product_dat_column_count(work_dirs):
    rows = _read_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if len(r) != 8]
    assert not bad, f"Rows with wrong column count (expected 8): {bad[:5]}"


def test_product_dat_codon_is_three_acgt(work_dirs):
    """Codons are three nucleotide characters; N is allowed for masked positions."""
    rows = _read_product(work_dirs)
    pattern = re.compile(r'^[ACGTN]{3}$')
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


def test_product_dat_product_count_non_negative(work_dirs):
    """product_count (col 4) is 0 for downstream-of-frameshift rows where the
    codon is disrupted; non-negative is the correct invariant."""
    rows = _read_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[3].lstrip('-').isdigit() or int(r[3]) < 0]
    assert not bad, f"Rows with product_count < 0 (col 4): {bad[:5]}"


def test_product_dat_amino_acid_single_char_or_stop(work_dirs):
    rows = _read_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if len(r[4]) != 1]
    assert not bad, f"Rows where amino_acid is not single char (col 5): {bad[:5]}"


def test_product_dat_pos_in_cds_positive(work_dirs):
    rows = _read_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if not r[5].isdigit() or int(r[5]) <= 0]
    assert not bad, f"Rows with pos_in_cds <= 0 (col 6): {bad[:5]}"


def test_product_dat_col7_pos_in_codon_1_2_3(work_dirs):
    rows = _read_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if r[6] not in ('1', '2', '3')]
    assert not bad, f"Rows with col 7 (dup pos_in_codon) not in {{1,2,3}}: {bad[:5]}"


def test_product_dat_downstream_of_frameshift_binary(work_dirs):
    rows = _read_product(work_dirs)
    bad = [i + 1 for i, r in enumerate(rows) if r[7] not in ('0', '1')]
    assert not bad, f"Rows with downstream_of_frameshift not 0/1 (col 8): {bad[:5]}"


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


def test_ann_vcf_has_some_cann_values(work_dirs):
    """CANN is only present on coding variants. Verify at least some records
    carry the CANN annotation (header presence verified separately)."""
    path = os.path.join(work_dirs['snpEff'], 'merged.ann.vcf.gz')
    values = bcftools_query_info(path, 'CANN')
    annotated = [v for v in values if v and v != '.']
    assert len(annotated) > 0, "No records have CANN in merged.ann.vcf.gz"


# ---------------------------------------------------------------------------
# Layer 2: Cross-output consistency
# ---------------------------------------------------------------------------

def test_strain_names_consistent_across_vcf_and_coverage(work_dirs):
    """VCF sample names and coverage.tsv strain columns must match exactly."""
    vcf_path = os.path.join(work_dirs['mergeVcfs'], 'merged.vcf.gz')
    vcf_strains = set(bcftools_samples(vcf_path))

    tsv_path = os.path.join(work_dirs['mergeCoverageBeds'], 'coverage.tsv')
    with open(tsv_path) as f:
        header = f.readline().rstrip('\n').split('\t')
    tsv_strains = set(header[3:])

    assert vcf_strains == tsv_strains, \
        f"VCF samples {vcf_strains} != coverage.tsv strains {tsv_strains}"


def test_strain_names_consistent_in_genomic_indel_db(work_dirs):
    """Strains in genomic_indels must be a subset of VCF sample names."""
    vcf_path = os.path.join(work_dirs['mergeVcfs'], 'merged.vcf.gz')
    vcf_strains = set(bcftools_samples(vcf_path))

    db_path = os.path.join(work_dirs['makeGenomicIndelDb'], 'genomicIndels.db')
    with sqlite3.connect(db_path) as conn:
        db_strains = {
            row[0]
            for row in conn.execute("SELECT DISTINCT strain FROM genomic_indels")
        }
    assert db_strains.issubset(vcf_strains), \
        f"genomic_indels strains not a subset of VCF samples: {db_strains - vcf_strains}"


def test_strain_names_consistent_in_coding_sequences(work_dirs):
    """VCF sample names must all be present in coding_sequences.db."""
    vcf_path = os.path.join(work_dirs['mergeVcfs'], 'merged.vcf.gz')
    vcf_strains = set(bcftools_samples(vcf_path))

    db_path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    with sqlite3.connect(db_path) as conn:
        db_strains = {
            row[0]
            for row in conn.execute("SELECT DISTINCT strain FROM coding_sequences")
        }
    assert vcf_strains.issubset(db_strains), \
        f"VCF strains not present in coding_sequences: {vcf_strains - db_strains}"


def test_variant_count_output_vcf_equals_ann_vcf(work_dirs):
    """snpEff annotates records but does not add or remove them."""
    output_vcf = os.path.join(work_dirs['processSeqVars'], 'output.vcf.gz')
    ann_vcf = os.path.join(work_dirs['snpEff'], 'merged.ann.vcf.gz')
    output_count = bcftools_record_count(output_vcf)
    ann_count = bcftools_record_count(ann_vcf)
    assert output_count == ann_count, \
        f"output.vcf.gz has {output_count} records but merged.ann.vcf.gz has {ann_count}"


def test_variation_feature_count_lte_vcf_record_count(work_dirs):
    """variationFeature.dat row count must be > 0 and <= VCF records x max transcripts per position."""
    vcf_path = os.path.join(work_dirs['processSeqVars'], 'output.vcf.gz')
    vcf_count = bcftools_record_count(vcf_path)

    vf_path = os.path.join(work_dirs['processSeqVars'], 'variationFeature.dat')
    with open(vf_path) as f:
        vf_count = sum(1 for _ in f)

    assert vf_count > 0, "variationFeature.dat is empty"
    assert vf_count <= vcf_count * 20, (
        f"variationFeature.dat rows ({vf_count}) suspiciously exceed "
        f"20x VCF records ({vcf_count})"
    )


def test_allele_dat_row_count_consistent_with_coding_variants(work_dirs):
    """allele.dat is only emitted for coding variants. Row count must be > 0
    and <= 4x coding variationFeature rows (conservatively: <= 4 distinct alleles per site)."""
    vf_rows = _read_variation_feature(work_dirs)
    coding_row_count = sum(1 for r in vf_rows if r[14] == '1')
    allele_row_count = len(_read_allele(work_dirs))
    assert allele_row_count > 0, "allele.dat is empty"
    assert allele_row_count <= coding_row_count * 4, (
        f"allele.dat rows ({allele_row_count}) far exceed "
        f"4x coding variationFeature rows ({coding_row_count})"
    )


def test_product_dat_transcripts_in_coding_sequences(work_dirs):
    """Every transcript_id in product.dat must exist in codingSequences.db."""
    db_path = os.path.join(work_dirs['makeCodingData'], 'codingSequences.db')
    with sqlite3.connect(db_path) as conn:
        known = {
            row[0]
            for row in conn.execute("SELECT DISTINCT transcript_id FROM coding_sequences")
        }

    prod_rows = _read_product(work_dirs)
    product_transcripts = {r[2] for r in prod_rows if r[2]}
    orphans = product_transcripts - known
    assert not orphans, \
        f"product.dat references transcripts not in codingSequences.db: {list(orphans)[:5]}"


def test_reference_strain_nonempty_and_consistent_with_vf(work_dirs):
    """variationFeature.dat col 4 (reference_strain) is the reference genome identifier
    used by GUS — not a VCF sample name. Verify it is non-empty and consistent across
    all rows (the cross-output check is that it never changes mid-file)."""
    vf_rows = _read_variation_feature(work_dirs)
    ref_strain = vf_rows[0][3]
    assert ref_strain, "reference_strain (col 4) is empty in first row of variationFeature.dat"
    inconsistent = [i + 1 for i, r in enumerate(vf_rows) if r[3] != ref_strain]
    assert not inconsistent, \
        f"reference_strain changes mid-file on rows: {inconsistent[:5]}"


def test_variation_feature_reference_strain_uniform(work_dirs):
    """All rows in variationFeature.dat must have the same reference_strain (col 4)."""
    vf_rows = _read_variation_feature(work_dirs)
    values = {r[3] for r in vf_rows}
    assert len(values) == 1, f"variationFeature.dat col 4 has multiple values: {values}"


def test_cann_present_in_some_output_vcf_records(work_dirs):
    """CANN is only set on coding variants. Verify at least some records
    have non-'.' CANN in output.vcf.gz (confirms coding annotation ran)."""
    path = os.path.join(work_dirs['processSeqVars'], 'output.vcf.gz')
    values = bcftools_query_info(path, 'CANN')
    annotated = [v for v in values if v and v != '.']
    assert len(annotated) > 0, \
        "No records have CANN in output.vcf.gz (expected at least some coding variants)"


def test_cann_and_ann_both_in_merged_ann_vcf(work_dirs):
    """merged.ann.vcf.gz must have ANN on all records and CANN on at least some."""
    path = os.path.join(work_dirs['snpEff'], 'merged.ann.vcf.gz')
    ann_values = bcftools_query_info(path, 'ANN')
    cann_values = bcftools_query_info(path, 'CANN')
    missing_ann = [i + 1 for i, v in enumerate(ann_values) if not v or v == '.']
    assert not missing_ann, f"Records missing ANN in merged.ann.vcf.gz: {missing_ann[:5]}"
    annotated_cann = [v for v in cann_values if v and v != '.']
    assert len(annotated_cann) > 0, \
        "No records have CANN in merged.ann.vcf.gz (expected at least some coding variants)"
