import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../bin'))

from makeConsensusFastaFromVcfAndBed import (
    load_coverage_bed,
    fill_gap,
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
    assert result == {'chr1': ([(0, 100), (200, 300)], [0, 200]), 'chr2': ([(50, 150)], [50])}


def test_load_coverage_bed_empty(tmp_path):
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    assert load_coverage_bed(str(bed)) == {}


# ── fill_gap ──────────────────────────────────────────────────────────────────

def test_fill_gap_all_covered():
    ref_seq = 'ACGTACGT'
    assert fill_gap(ref_seq, 0, 8, [(0, 8)], [0]) == 'ACGTACGT'


def test_fill_gap_all_uncovered():
    assert fill_gap('ACGTACGT', 0, 8, [], []) == 'NNNNNNNN'


def test_fill_gap_partial_coverage():
    # covered [2,5), uncovered [0,2) and [5,8)
    ref_seq = 'ACGTACGT'
    result = fill_gap(ref_seq, 0, 8, [(2, 5)], [2])
    assert result == 'NN' + 'GTA' + 'NNN'


def test_fill_gap_subrange_with_partial_overlap():
    # fill [1,6), coverage at [2,5)
    ref_seq = 'ACGTACGT'
    result = fill_gap(ref_seq, 1, 6, [(2, 5)], [2])
    assert result == 'N' + 'GTA' + 'N'


def test_fill_gap_empty_range():
    assert fill_gap('ACGT', 3, 3, [(0, 4)], [0]) == ''


def test_fill_gap_multiple_intervals():
    # covered [0,3) and [5,8), uncovered [3,5)
    ref_seq = 'ACGTACGT'
    result = fill_gap(ref_seq, 0, 8, [(0, 3), (5, 8)], [0, 5])
    assert result == 'ACG' + 'NN' + 'CGT'


# ── build_consensus ───────────────────────────────────────────────────────────

def test_build_consensus_all_covered_no_variants():
    ref_seq = 'ACGTACGT'
    vcf = FakeVcf({})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(0, 8)], [0])
    assert seq == 'ACGTACGT'


def test_build_consensus_all_uncovered_no_variants():
    ref_seq = 'ACGTACGT'
    vcf = FakeVcf({})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [], [])
    assert seq == 'NNNNNNNN'


def test_build_consensus_snp_covered():
    # SNP at pos 5 (1-based), A->T hom, fully covered
    ref_seq = 'ACGTACGT'
    record = FakeVcfRecord('chr1', 5, 'A', ['T'], 'T/T')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(0, 8)], [0])
    assert seq == 'ACGT' + 'T' + 'CGT'


def test_build_consensus_snp_het_iupac():
    # het A/T at pos 1 → IUPAC W
    ref_seq = 'ACGT'
    record = FakeVcfRecord('chr1', 1, 'A', ['T'], 'A/T')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 4, ref_seq, vcf, [(0, 4)], [0])
    assert seq[0] == 'W'
    assert seq[1:] == 'CGT'


def test_build_consensus_snp_applied_even_when_bed_uncovered():
    # Variant positions are authoritative: FreeBayes enforces --min-coverage, so
    # a called SNP is applied regardless of BED coverage. The empty BED masks
    # only the non-variant gap after it. See build_consensus docstring.
    ref_seq = 'ACGT'
    record = FakeVcfRecord('chr1', 1, 'A', ['T'], 'T/T')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 4, ref_seq, vcf, [], [])
    assert seq == 'TNNN'


def test_build_consensus_hom_insertion_covered():
    # Hom insertion at pos 3: REF=G, ALT=GCC (len 3 vs 1)
    ref_seq = 'ACGTACGT'
    record = FakeVcfRecord('chr1', 3, 'G', ['GCC'], 'GCC/GCC')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(0, 8)], [0])
    assert seq == 'AC' + 'GCC' + 'TACGT'


def test_build_consensus_hom_deletion_covered():
    # Hom deletion at pos 1: REF=ACG, ALT=A (del of CG)
    ref_seq = 'ACGTACGT'
    record = FakeVcfRecord('chr1', 1, 'ACG', ['A'], 'A/A')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(0, 8)], [0])
    assert seq == 'A' + 'TACGT'


def test_build_consensus_het_indel_01_emits_ref():
    # Het 0/1: one allele is REF → emit REF
    ref_seq = 'ACGTACGT'
    record = FakeVcfRecord('chr1', 2, 'CG', ['C'], 'CG/C')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(0, 8)], [0])
    assert seq == 'A' + 'CG' + 'TACGT'


def test_build_consensus_het_indel_12_emits_x_times_ref_len():
    # Het 1/2: both non-ref, REF=ACG (len 3) → XXX
    ref_seq = 'ACGTACGT'
    record = FakeVcfRecord('chr1', 1, 'ACG', ['AC', 'A'], 'AC/A')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(0, 8)], [0])
    assert seq == 'XXX' + 'TACGT'


def test_build_consensus_gap_before_variant_partially_covered():
    # variant at pos 5 (1-based=0-based 4), coverage only [2,8)
    # gap [0,4): [0,2) uncovered → NN, [2,4) covered → ref GT
    ref_seq = 'ACGTACGT'
    record = FakeVcfRecord('chr1', 5, 'A', ['T'], 'T/T')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(2, 8)], [2])
    assert seq == 'NN' + 'GT' + 'T' + 'CGT'


def test_build_consensus_gap_after_last_variant():
    # variant at pos 1, coverage [0,8), fill [1,8) after it
    ref_seq = 'ACGTACGT'
    record = FakeVcfRecord('chr1', 1, 'A', ['T'], 'T/T')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(0, 8)], [0])
    assert seq == 'T' + 'CGTACGT'


def test_build_consensus_dot_gt_gives_n():
    ref_seq = 'ACGT'
    record = FakeVcfRecord('chr1', 1, 'A', ['T'], './.')
    vcf = FakeVcf({'chr1': [record]})
    seq = build_consensus('chr1', 4, ref_seq, vcf, [(0, 4)], [0])
    assert seq == 'N' + 'CGT'


def test_build_consensus_overlapping_records_skipped():
    # Second record at pos 2 (0-based 1) overlaps with first record REF=ACG (spans 0-2).
    # The overlapping record should be skipped; output length must equal chrom_len.
    ref_seq = 'ACGTACGT'
    r1 = FakeVcfRecord('chr1', 1, 'ACG', ['TTT'], 'TTT/TTT')  # hom, pos 0-2
    r2 = FakeVcfRecord('chr1', 2, 'C', ['G'], 'G/G')           # overlaps — pos 1, inside r1
    vcf = FakeVcf({'chr1': [r1, r2]})
    seq = build_consensus('chr1', 8, ref_seq, vcf, [(0, 8)], [0])
    assert len(seq) == 8
    assert seq == 'TTT' + 'TACGT'


# ── write_fasta ───────────────────────────────────────────────────────────────

def test_write_fasta(tmp_path):
    out = tmp_path / "out.fa"
    with open(out, 'w') as fh:
        write_fasta(fh, 'chr1', 'ACGT' * 20)
    lines = out.read_text().splitlines()
    assert lines[0] == '>chr1'
    assert ''.join(lines[1:]) == 'ACGT' * 20


def test_write_fasta_sample_prefixed_defline(tmp_path):
    # main() composes the defline as "<sample>_<sequenceId>"; downstream
    # makeCodingData.jl strips the "<strain>_" prefix to recover the bare
    # sequenceId. Lock the convention here.
    out = tmp_path / "out.fa"
    sample, chrom = 'strainA', 'chr1'
    with open(out, 'w') as fh:
        write_fasta(fh, f'{sample}_{chrom}', 'ACGT')
    assert out.read_text().splitlines()[0] == '>strainA_chr1'
