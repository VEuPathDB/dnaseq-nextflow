import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../bin'))

from mergeVariantsByLocation import _recombine, _fix_single, _flush


def _make_row(chrom, pos, ref, alt, fmt, sample):
    return [chrom, str(pos), '.', ref, alt, '.', 'PASS', '.', fmt, sample]


# ── _fix_single ───────────────────────────────────────────────────────────────

def test_fix_single_no_star_passthrough():
    row = _make_row('chr1', 100, 'A', 'C', 'GT:AD', '0/1:10,5')
    assert _fix_single(row) == row


def test_fix_single_star_only_dropped():
    row = _make_row('chr1', 100, 'A', '*', 'GT:AD', '0/1:10,.')
    assert _fix_single(row) is None


def test_fix_single_real_alt_plus_star_remaps_gt():
    # ALT=C,* → GT 0/1 stays 0/1; * slot stripped from AD
    row = _make_row('chr1', 100, 'A', 'C,*', 'GT:AD', '0/1:10,5,.')
    result = _fix_single(row)
    assert result is not None
    assert result[4] == 'C'
    sample = result[9].split(':')
    assert sample[0] == '0/1'   # GT unchanged (C is still allele 1)
    assert sample[1] == '10,5'  # AD: REF + C; * slot removed


# ── _recombine: GT ────────────────────────────────────────────────────────────

def test_recombine_gt_1_2():
    # Standard het multi-allelic: GT 1/2 split across two rows
    row0 = _make_row('chr1', 83041, 'A', 'C,*', 'GT', '1/2')
    row1 = _make_row('chr1', 83041, 'A', '*,G', 'GT', '1/2')
    result = _recombine([row0, row1])
    assert result is not None
    assert result[4] == 'C,G'
    assert result[9] == '1/2'


def test_recombine_gt_0_1_and_0_2():
    # One row hom-ref for its real allele, other is alt
    row0 = _make_row('chr1', 200, 'T', 'A,*', 'GT', '0/1')
    row1 = _make_row('chr1', 200, 'T', '*,G', 'GT', '0/2')
    result = _recombine([row0, row1])
    assert result is not None
    assert result[4] == 'A,G'
    gt_alleles = set(result[9].split('/'))
    assert '0' in gt_alleles


# ── _recombine: AD (Number=R) ─────────────────────────────────────────────────

def test_recombine_ad_merges_real_allele_depths():
    # Issue #9: row0 AD=0,11,. and row1 AD=0,.,12 → combined 0,11,12
    row0 = _make_row('chr1', 83041, 'A', 'C,*', 'GT:AD', '1/2:0,11,.')
    row1 = _make_row('chr1', 83041, 'A', '*,G', 'GT:AD', '1/2:0,.,12')
    result = _recombine([row0, row1])
    assert result is not None
    sample_fields = result[9].split(':')
    assert sample_fields[1] == '0,11,12'


def test_recombine_ad_star_first_row():
    # * is the first allele in row1: ALT=*,G
    row0 = _make_row('chr1', 300, 'C', 'T,*', 'GT:AD', '1/2:0,8,.')
    row1 = _make_row('chr1', 300, 'C', '*,A', 'GT:AD', '1/2:0,.,5')
    result = _recombine([row0, row1])
    assert result is not None
    sample_fields = result[9].split(':')
    assert sample_fields[1] == '0,8,5'


# ── _recombine: AO (Number=A) ─────────────────────────────────────────────────

def test_recombine_ao_merges_real_allele_counts():
    row0 = _make_row('chr1', 83041, 'A', 'C,*', 'GT:AO', '1/2:11,.')
    row1 = _make_row('chr1', 83041, 'A', '*,G', 'GT:AO', '1/2:.,12')
    result = _recombine([row0, row1])
    assert result is not None
    sample_fields = result[9].split(':')
    assert sample_fields[1] == '11,12'


def test_recombine_ao_star_first():
    row0 = _make_row('chr1', 400, 'G', 'T,*', 'GT:AO', '1/2:7,.')
    row1 = _make_row('chr1', 400, 'G', '*,C', 'GT:AO', '1/2:.,3')
    result = _recombine([row0, row1])
    assert result is not None
    sample_fields = result[9].split(':')
    assert sample_fields[1] == '7,3'


# ── _recombine: Number=1 fields unchanged ─────────────────────────────────────

def test_recombine_dp_kept_from_row0():
    row0 = _make_row('chr1', 83041, 'A', 'C,*', 'GT:AD:DP', '1/2:0,11,.:30')
    row1 = _make_row('chr1', 83041, 'A', '*,G', 'GT:AD:DP', '1/2:0,.,12:25')
    result = _recombine([row0, row1])
    assert result is not None
    sample_fields = result[9].split(':')
    assert sample_fields[2] == '30'  # DP from row 0


# ── _recombine: combined AD + AO + DP ─────────────────────────────────────────

def test_recombine_combined_fields():
    row0 = _make_row('chr1', 83041, 'A', 'C,*', 'GT:AD:AO:DP', '1/2:0,11,.:11,.:30')
    row1 = _make_row('chr1', 83041, 'A', '*,G', 'GT:AD:AO:DP', '1/2:0,.,12:.,12:25')
    result = _recombine([row0, row1])
    assert result is not None
    assert result[4] == 'C,G'
    sample_fields = result[9].split(':')
    assert sample_fields[0] == '1/2'    # GT
    assert sample_fields[1] == '0,11,12'  # AD (Number=R)
    assert sample_fields[2] == '11,12'    # AO (Number=A)
    assert sample_fields[3] == '30'       # DP (Number=1)


# ── _flush integration ────────────────────────────────────────────────────────

def test_flush_single_row_no_star_passthrough():
    row = '\t'.join(_make_row('chr1', 100, 'A', 'C', 'GT:AD', '0/1:10,5')) + '\n'
    out = _flush([row])
    assert len(out) == 1
    assert 'C' in out[0]


def test_flush_two_rows_recombined():
    row0 = '\t'.join(_make_row('chr1', 83041, 'A', 'C,*', 'GT:AD', '1/2:0,11,.')) + '\n'
    row1 = '\t'.join(_make_row('chr1', 83041, 'A', '*,G', 'GT:AD', '1/2:0,.,12')) + '\n'
    out = _flush([row0, row1])
    assert len(out) == 1
    fields = out[0].rstrip('\n').split('\t')
    assert fields[4] == 'C,G'
    assert fields[9].split(':')[1] == '0,11,12'
