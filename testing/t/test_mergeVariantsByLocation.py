import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../bin'))

from mergeVariantsByLocation import _recombine, _fix_single, _flush

# Real data format from FreeBayes/filterAndSplitVcf pipeline.
# All rows below are derived from LmjF.01 (chr1) runs in the WhitePaper workflow.
FMT = 'GT:DP:AD:RO:QR:AO:QA:GL'


def _make_row(chrom, pos, ref, alt, fmt, sample):
    return [chrom, str(pos), '.', ref, alt, '.', 'PASS', '.', fmt, sample]


# ── _fix_single ───────────────────────────────────────────────────────────────

def test_fix_single_no_star_passthrough():
    # Real data: LmjF.01:420 simple hom-alt SNP, no * in ALT
    row = _make_row('LmjF.01', 420, 'C', 'G', FMT,
                    '1/1:18:0,18:0:0:18:1257:-97.0233,-5.41854,0')
    assert _fix_single(row) == row


def test_fix_single_star_only_dropped():
    # * as sole ALT: row is dropped entirely.
    # Does not appear in real data (bcftools norm -a always preserves a real alt),
    # but the code path exists and must not silently pass through.
    row = _make_row('LmjF.01', 100, 'A', '*', 'GT:AD', '0/1:10,.')
    assert _fix_single(row) is None


def test_fix_single_gt_0_2_becomes_0_1():
    # Real data: LmjF.01:58460 C→G,* with GT=0/2 (Case B/C — solo row after
    # bcftools norm -a splits a complex call; partner has different REF).
    # Allele 2 points to *; strips * and remaps 2→1 throughout.
    row = _make_row('LmjF.01', 58460, 'C', 'G,*', FMT,
                    '0/2:79:37,31,.:37:2510:31,.:2101,.:-139.304,0,-174.79,.,.,.')
    result = _fix_single(row)
    assert result is not None
    assert result[4] == 'G'
    parts = result[9].split(':')
    assert parts[0] == '0/1'    # GT: 0→0, 2→1
    assert parts[2] == '37,31'  # AD (Number=R): * slot removed
    assert parts[5] == '31'     # AO (Number=A): * slot removed


def test_fix_single_gt_1_2_collapses_to_1_1():
    # Real data: LmjF.01:68344 GTA→G,* with GT=1/2.
    # Both allele 1 (real alt) and allele 2 (*) remap to new index 1 → GT=1/1.
    # Arises when bcftools norm -a partner (adjacent position) was filtered out.
    row = _make_row('LmjF.01', 68344, 'GTA', 'G,*', FMT,
                    '1/2:15:0,9,.:0:0:9,.:572,.:-71.5214,-29.4946,-26.7853,.,.,.')
    result = _fix_single(row)
    assert result is not None
    assert result[4] == 'G'
    parts = result[9].split(':')
    assert parts[0] == '1/1'   # GT: 1→1, 2→1
    assert parts[2] == '0,9'   # AD (Number=R): * slot removed


def test_fix_single_gt_1_1_unchanged():
    # Real data: LmjF.01:141874 G→C,* with GT=1/1 (hom-alt, * not in GT).
    # GT stays 1/1; * is stripped from AD/AO.
    row = _make_row('LmjF.01', 141874, 'G', 'C,*', FMT,
                    '1/1:264:42,202,.:42:2847:120,.:8135,.:-751.913,-171.602,-339.113,.,.,.')
    result = _fix_single(row)
    assert result is not None
    assert result[4] == 'C'
    parts = result[9].split(':')
    assert parts[0] == '1/1'
    assert parts[2] == '42,202'  # AD (Number=R): * slot removed
    assert parts[5] == '120'     # AO (Number=A): * slot removed


# ── _recombine ────────────────────────────────────────────────────────────────

def test_recombine_gt_1_2_and_2_1():
    # Real data: LmjF.01:8962 T→{TA,TAA} split by bcftools norm -a.
    # Row0 carries TA (GT allele 1 points to it), row1 carries TAA (GT allele 1
    # points to it after * is at index 2). Combined GT=1/2, AD and AO merged.
    row0 = _make_row('LmjF.01', 8962, 'T', 'TA,*', FMT,
                     '1/2:12:0,7,.:0:0:7,.:406,.:-55.505,-22.5762,-20.469,.,.,.')
    row1 = _make_row('LmjF.01', 8962, 'T', 'TAA,*', FMT,
                     '2/1:12:0,5,.:0:0:5,.:272,.:-55.505,-33.4657,-31.9606,.,.,.')
    result = _recombine([row0, row1])
    assert result is not None
    assert result[4] == 'TA,TAA'
    parts = result[9].split(':')
    assert parts[0] == '1/2'
    assert parts[2] == '0,7,5'  # AD (Number=R): REF + TA + TAA
    assert parts[5] == '7,5'    # AO (Number=A): TA + TAA
    assert parts[1] == '12'     # DP (Number=1): kept from row0


def test_recombine_het_snp_pair():
    # Real data: LmjF.01:141701 T→{C,G} — two SNP alts split by bcftools norm -a.
    # Confirms recombination works for SNPs, not just indels.
    row0 = _make_row('LmjF.01', 141701, 'T', 'C,*', FMT,
                     '1/2:56:13,26,.:13:863:26,.:1788,.:-138.029,-13.1137,-69.6537,.,.,.')
    row1 = _make_row('LmjF.01', 141701, 'T', 'G,*', FMT,
                     '2/1:56:13,17,.:13:863:17,.:1185,.:-138.029,-64.4548,-123.685,.,.,.')
    result = _recombine([row0, row1])
    assert result is not None
    assert result[4] == 'C,G'
    parts = result[9].split(':')
    assert parts[0] == '1/2'
    assert parts[2] == '13,26,17'  # AD (Number=R): REF + C + G
    assert parts[5] == '26,17'    # AO (Number=A): C + G
    assert parts[1] == '56'       # DP (Number=1): kept from row0


def test_recombine_gt_2_2_fallback_positional():
    # Real data: LmjF.01:66728 A→{AGT,T} where both rows have GT=2/2.
    # Allele 2 in each row points to * — neither row's GT references the real alt.
    # Falls back to positional assignment: row0→allele 1, row1→allele 2 → GT=1/2.
    row0 = _make_row('LmjF.01', 66728, 'A', 'AGT,*', FMT,
                     '2/2:21:0,17,.:0:0:17,.:1014,.:-86.8211,-5.11751,0,.,.,.')
    row1 = _make_row('LmjF.01', 66728, 'A', 'T,*', FMT,
                     '2/2:21:0,17,.:0:0:17,.:1014,.:-86.8211,-5.11751,0,.,.,.')
    result = _recombine([row0, row1])
    assert result is not None
    assert result[4] == 'AGT,T'
    parts = result[9].split(':')
    assert parts[0] == '1/2'       # fallback: row0→1, row1→2
    assert parts[2] == '0,17,17'   # AD (Number=R): REF + AGT + T
    assert parts[5] == '17,17'     # AO (Number=A): AGT + T


# ── _flush integration ────────────────────────────────────────────────────────

def test_flush_single_no_star_passthrough():
    # Single row with no * passes through byte-for-byte unchanged
    row = '\t'.join(_make_row('LmjF.01', 420, 'C', 'G', FMT,
                              '1/1:18:0,18:0:0:18:1257:-97.0233,-5.41854,0')) + '\n'
    out = _flush([row])
    assert len(out) == 1
    assert out[0] == row


def test_flush_real_recombine_pair():
    # Real data pair LmjF.01:8962 — two * rows recombined into one
    row0 = '\t'.join(_make_row('LmjF.01', 8962, 'T', 'TA,*', FMT,
                               '1/2:12:0,7,.:0:0:7,.:406,.:-55.505,-22.5762,-20.469,.,.,.'))\
           + '\n'
    row1 = '\t'.join(_make_row('LmjF.01', 8962, 'T', 'TAA,*', FMT,
                               '2/1:12:0,5,.:0:0:5,.:272,.:-55.505,-33.4657,-31.9606,.,.,.'))\
           + '\n'
    out = _flush([row0, row1])
    assert len(out) == 1
    fields = out[0].rstrip('\n').split('\t')
    assert fields[4] == 'TA,TAA'
    parts = fields[9].split(':')
    assert parts[0] == '1/2'
    assert parts[2] == '0,7,5'


def test_flush_real_fix_single_gt_remapped():
    # Single * row (Case B/C) — GT=0/2 remapped to GT=0/1, * stripped
    row = '\t'.join(_make_row('LmjF.01', 58460, 'C', 'G,*', FMT,
                              '0/2:79:37,31,.:37:2510:31,.:2101,.:-139.304,0,-174.79,.,.,.'))\
          + '\n'
    out = _flush([row])
    assert len(out) == 1
    fields = out[0].rstrip('\n').split('\t')
    assert fields[4] == 'G'
    assert fields[9].split(':')[0] == '0/1'
