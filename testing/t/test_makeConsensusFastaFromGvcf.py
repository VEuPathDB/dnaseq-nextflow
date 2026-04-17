#!/usr/bin/env python3
import sys, os, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../bin'))

from makeConsensusFastaFromGvcf import build_consensus

class FakeVariant:
    def __init__(self, pos, ref, alts, dp, gt_bases_val):
        self.POS = pos
        self.REF = ref
        self.ALT = alts
        self._dp = dp
        self._gt_bases = gt_bases_val
        self.FORMAT = ['GT', 'DP']

    def format(self, key):
        if key == 'DP':
            return [[self._dp]]
        return None

    @property
    def gt_bases(self):
        return [self._gt_bases]

class FakeVCF:
    def __init__(self, records):
        self._records = records
    def __call__(self, chrom):
        return iter(self._records)

class TestBuildConsensus(unittest.TestCase):

    def _run(self, chrom_len, ref_seq, records, min_coverage=5):
        vcf = FakeVCF(records)
        return build_consensus('chr1', chrom_len, ref_seq, vcf, min_coverage)

    def test_hom_deletion_emits_ref_bases(self):
        # Hom deletion: REF=ATGC at pos 3 (1-based), ALT=A
        ref_seq = 'XXATGCXX'
        rec = FakeVariant(pos=3, ref='ATGC', alts=['A'], dp=20, gt_bases_val='A/A')
        result = self._run(len(ref_seq), ref_seq, [rec])
        self.assertEqual(len(result), len(ref_seq), 'output must be reference-length')
        self.assertEqual(result[2:6], 'ATGC', 'hom deletion emits ref bases')

    def test_het_deletion_emits_ref_bases(self):
        # Het deletion: REF=ATGC at pos 3 (1-based), ALT=A
        ref_seq = 'XXATGCXX'
        rec = FakeVariant(pos=3, ref='ATGC', alts=['A'], dp=20, gt_bases_val='ATGC/A')
        result = self._run(len(ref_seq), ref_seq, [rec])
        self.assertEqual(len(result), len(ref_seq), 'output must be reference-length')
        self.assertEqual(result[2:6], 'ATGC', 'het deletion emits ref bases')

    def test_hom_insertion_emits_ref_bases(self):
        # Hom insertion: REF=A at pos 3 (1-based), ALT=ATGC
        ref_seq = 'XXAXX'
        rec = FakeVariant(pos=3, ref='A', alts=['ATGC'], dp=20, gt_bases_val='ATGC/ATGC')
        result = self._run(len(ref_seq), ref_seq, [rec])
        self.assertEqual(len(result), len(ref_seq), 'output must be reference-length')
        self.assertEqual(result[2], 'A', 'hom insertion emits ref base')

    def test_low_coverage_emits_n(self):
        ref_seq = 'XXATGCXX'
        rec = FakeVariant(pos=3, ref='ATGC', alts=['A'], dp=2, gt_bases_val='A/A')
        result = self._run(len(ref_seq), ref_seq, [rec], min_coverage=5)
        self.assertEqual(result[2:6], 'NNNN', 'low coverage emits N')

if __name__ == '__main__':
    unittest.main()
