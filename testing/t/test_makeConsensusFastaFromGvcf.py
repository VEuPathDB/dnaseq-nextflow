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
