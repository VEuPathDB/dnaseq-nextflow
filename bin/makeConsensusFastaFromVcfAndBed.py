#!/usr/bin/env python3
"""
Build a consensus FASTA from a standard VCF (FreeBayes / bcftools) plus a
BED file that defines covered regions.

Coverage at VCF record positions is not re-checked: FreeBayes already enforces
--min-coverage, so every called variant has sufficient depth.  The BED is used
only to fill non-variant gaps between records with ref sequence or N.
"""

import argparse
import bisect
import subprocess
from collections import defaultdict

# IUPAC ambiguity codes keyed by frozenset of bases
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


def load_coverage_bed(bed_path):
    """
    Parse a BED file (plain or gzipped) and return a dict mapping
    chrom -> (intervals, starts) where intervals is a sorted list of
    (start, end) tuples (0-based, half-open) and starts is the precomputed
    list of interval start coordinates for binary search.
    """
    import gzip as _gzip
    opener = _gzip.open if bed_path.endswith('.gz') else open
    raw = defaultdict(list)
    with opener(bed_path, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            raw[chrom].append((start, end))
    result = {}
    for chrom, ivs in raw.items():
        ivs.sort()
        result[chrom] = (ivs, [iv[0] for iv in ivs])
    return result


def fill_gap(ref_seq, gap_start, gap_end, intervals, starts):
    """Fill [gap_start, gap_end) with ref bases where covered, N elsewhere. No per-base loop."""
    if gap_start >= gap_end:
        return ''
    if not intervals:
        return 'N' * (gap_end - gap_start)

    segments = []
    pos = gap_start
    idx = max(bisect.bisect_right(starts, gap_start) - 1, 0)

    for iv_start, iv_end in intervals[idx:]:
        if iv_start >= gap_end:
            break
        fill_to = min(iv_start, gap_end)
        if pos < fill_to:
            segments.append('N' * (fill_to - pos))
            pos = fill_to
        cov_start = max(iv_start, pos)
        cov_end = min(iv_end, gap_end)
        if cov_start < cov_end:
            segments.append(ref_seq[cov_start:cov_end])
            pos = cov_end

    if pos < gap_end:
        segments.append('N' * (gap_end - pos))
    return ''.join(segments)


def build_consensus(chrom_name, chrom_len, ref_seq, vcf, intervals, starts):
    """
    Walk VCF records for one chromosome and assemble the consensus sequence.

    Parameters
    ----------
    chrom_name : str
    chrom_len  : int
    ref_seq    : str   – full reference sequence for the chromosome (0-based)
    vcf        : callable – vcf(chrom_name) returns an iterator of VCF records
    intervals  : list of (start, end) tuples (0-based, half-open) from BED, sorted by start
    starts     : list of interval start coordinates (precomputed for binary search)

    Logic
    -----
    SNPs        → IUPAC code derived from GT alleles
    Hom indel   → emit the called allele sequence
    Het 0/1     → emit REF (one allele is reference)
    Het 1/2     → 'X' × len(REF) (both alleles non-ref, ambiguous)
    Missing GT  → 'N' × len(REF)
    Gaps        → fill_gap (ref or N based on BED coverage)

    Coverage at VCF positions is not re-checked against the BED: FreeBayes
    already enforces --min-coverage, so every record in the VCF has sufficient
    depth.  The BED is used only to fill non-variant gaps with ref or N.
    """
    segments = []
    ref_pos = 0   # current 0-based reference coordinate

    for v in vcf(chrom_name):
        pos = v.POS - 1   # VCF 1-based → 0-based
        if pos < ref_pos:   # overlapping record — skip
            continue

        # Fill any gap between the last processed position and this record
        if pos > ref_pos:
            segments.append(fill_gap(ref_seq, ref_pos, pos, intervals, starts))
            ref_pos = pos

        # GT string for the single sample
        gt_str = v.gt_bases[0] if v.gt_bases else './.'

        if '.' in gt_str:
            segments.append('N' * len(v.REF))
            ref_pos = pos + len(v.REF)
            continue

        alleles = list(dict.fromkeys(gt_str.replace('|', '/').split('/')))  # unique, ordered

        is_indel = len(v.REF) > 1 or any(len(a) > 1 for a in alleles)

        if not is_indel and all(len(a) == 1 for a in alleles):
            # ── SNP (or hom-ref) ──────────────────────────────────────────
            base = IUPAC.get(frozenset(alleles), 'N')
            segments.append(base)
            ref_pos = pos + 1

        elif len(alleles) == 1:
            # ── Homozygous indel ─────────────────────────────────────────
            segments.append(alleles[0])
            ref_pos = pos + len(v.REF)

        else:
            # ── Heterozygous indel ────────────────────────────────────────
            if v.REF in alleles:
                # 0/1 – one allele is REF; emit REF unchanged
                segments.append(v.REF)
            else:
                # 1/2 – both alleles non-ref
                ref_len = len(v.REF)
                if all(len(a) == ref_len for a in alleles):
                    # MNP: same-length alleles — decompose per position into IUPAC
                    out = []
                    for i in range(ref_len):
                        bases = frozenset(a[i] for a in alleles)
                        out.append(IUPAC.get(bases, 'N'))
                    segments.append(''.join(out))
                else:
                    # True indel with length change; mask with X
                    segments.append('X' * ref_len)
            ref_pos = pos + len(v.REF)

    # Fill remaining reference positions after the last VCF record
    if ref_pos < chrom_len:
        segments.append(fill_gap(ref_seq, ref_pos, chrom_len, intervals, starts))

    return ''.join(segments)


def write_fasta(out_fh, name, seq, line_len=60):
    out_fh.write(f'>{name}\n')
    for i in range(0, len(seq), line_len):
        out_fh.write(seq[i:i + line_len] + '\n')


def get_chrom_seq(ref_fasta, chrom_name):
    """Load a single chromosome sequence from an indexed FASTA."""
    result = subprocess.run(
        ['samtools', 'faidx', ref_fasta, chrom_name],
        capture_output=True, text=True, check=True,
    )
    return ''.join(result.stdout.splitlines()[1:]).upper()


def main():
    parser = argparse.ArgumentParser(
        description='Build a consensus FASTA from a VCF + BED coverage file.')
    parser.add_argument('-v', '--vcf',    required=True,
                        help='VCF or VCF.gz with tabix index (.tbi)')
    parser.add_argument('-r', '--ref',   required=True,
                        help='Reference FASTA (indexed with samtools faidx)')
    parser.add_argument('-f', '--fai',   required=True,
                        help='Reference FASTA .fai index')
    parser.add_argument('-b', '--bed',   required=True,
                        help='BED file of covered regions (0-based, half-open)')
    parser.add_argument('-o', '--output', required=True)
    args = parser.parse_args()

    chroms = []
    with open(args.fai) as fh:
        for line in fh:
            parts = line.split('\t')
            chroms.append((parts[0], int(parts[1])))

    from cyvcf2 import VCF  # deferred so tests can import without cyvcf2 installed

    coverage = load_coverage_bed(args.bed)
    vcf = VCF(args.vcf)

    with open(args.output, 'w') as out:
        for chrom_name, chrom_len in chroms:
            ref_seq = get_chrom_seq(args.ref, chrom_name)
            intervals, starts = coverage.get(chrom_name, ([], []))
            seq = build_consensus(chrom_name, chrom_len, ref_seq, vcf, intervals, starts)
            write_fasta(out, chrom_name, seq)


if __name__ == '__main__':
    main()
