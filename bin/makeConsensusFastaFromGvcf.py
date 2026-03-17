#!/usr/bin/env python3
import argparse
import subprocess
from cyvcf2 import VCF

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


def get_chrom_seq(ref_fasta, chrom_name):
    """Load a single chromosome sequence from an indexed FASTA."""
    result = subprocess.run(
        ['samtools', 'faidx', ref_fasta, chrom_name],
        capture_output=True, text=True, check=True,
    )
    return ''.join(result.stdout.splitlines()[1:]).upper()


def build_consensus(chrom_name, chrom_len, ref_seq, vcf, min_coverage):
    """
    Walk through g.vcf records for one chromosome and assemble the consensus.

    REF blocks   → reference bases from ref_seq (or N if dp < min_coverage)
    SNPs         → IUPAC code derived from GT alleles
    Indels       → homozygous: the called allele sequence (length may differ
                   from REF, so the output FASTA diverges from the reference);
                   heterozygous: N × len(REF) since two different-length
                   alleles cannot be collapsed into one sequence
    Gaps in VCF  → N (no coverage)
    """
    segments = []
    ref_pos = 0   # current position in reference coordinates (0-based)

    for v in vcf(chrom_name):
        pos = v.POS - 1   # VCF 1-based → 0-based

        # Fill any gap before this record with Ns (positions absent from g.vcf)
        if pos > ref_pos:
            segments.append('N' * (pos - ref_pos))
            ref_pos = pos

        dp_arr = v.format('DP')
        dp = int(dp_arr[0][0]) if dp_arr is not None else 0

        # ── REF block ──────────────────────────────────────────────────────
        # ALT is absent, '.', or a symbolic allele (<*>, <NON_REF>, …)
        if not v.ALT or all(a == '.' or a.startswith('<') for a in v.ALT):
            end = (v.INFO.get('END') or v.POS) - 1   # 0-based inclusive
            span = end - pos + 1
            if dp >= min_coverage:
                segments.append(ref_seq[pos:end + 1])
            else:
                segments.append('N' * span)
            ref_pos = end + 1
            continue

        # ── Variant record ─────────────────────────────────────────────────
        if dp < min_coverage:
            # Mask the REF span and advance
            segments.append('N' * len(v.REF))
            ref_pos = pos + len(v.REF)
            continue

        if 'GT' not in v.FORMAT:
            segments.append('N' * len(v.REF))
            ref_pos = pos + len(v.REF)
            continue

        gt_str = v.gt_bases[0]   # e.g. 'A/G' or 'A|G' for the single sample
        if '.' in gt_str:
            segments.append('N' * len(v.REF))
            ref_pos = pos + len(v.REF)
            continue

        alleles = list(dict.fromkeys(gt_str.replace('|', '/').split('/')))   # unique, ordered

        if all(len(a) == 1 for a in alleles):
            # ── SNP (or hom-ref call) ──
            base = IUPAC.get(frozenset(alleles), 'N')
            segments.append(base)
            ref_pos = pos + 1

        elif len(alleles) == 1:
            # ── Homozygous indel: emit the actual allele sequence ──
            segments.append(alleles[0])
            ref_pos = pos + len(v.REF)   # advance past REF span in reference coords

        else:
            # ── Heterozygous indel: two different-length alleles, mask ──
            segments.append('N' * len(v.REF))
            ref_pos = pos + len(v.REF)

    # Fill any remaining reference positions that had no coverage
    if ref_pos < chrom_len:
        segments.append('N' * (chrom_len - ref_pos))

    return ''.join(segments)


def write_fasta(out_fh, name, seq, line_len=60):
    out_fh.write(f'>{name}\n')
    for i in range(0, len(seq), line_len):
        out_fh.write(seq[i:i + line_len] + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Build a consensus FASTA from a g.vcf, applying IUPAC '
                    'codes for SNPs and actual allele sequences for indels.')
    parser.add_argument('-g', '--gvcf',         required=True,
                        help='g.vcf.gz with tabix index (.tbi)')
    parser.add_argument('-r', '--ref',          required=True,
                        help='Reference FASTA (indexed with samtools faidx)')
    parser.add_argument('-f', '--fai',          required=True,
                        help='Reference FASTA .fai index')
    parser.add_argument('-mc', '--min-coverage', required=True, type=int,
                        dest='min_coverage')
    parser.add_argument('-o', '--output',       required=True)
    args = parser.parse_args()

    chroms = []
    with open(args.fai) as fh:
        for line in fh:
            parts = line.split('\t')
            chroms.append((parts[0], int(parts[1])))

    vcf = VCF(args.gvcf)

    with open(args.output, 'w') as out:
        for chrom_name, chrom_len in chroms:
            ref_seq = get_chrom_seq(args.ref, chrom_name)
            seq = build_consensus(chrom_name, chrom_len, ref_seq, vcf, args.min_coverage)
            write_fasta(out, chrom_name, seq)


if __name__ == '__main__':
    main()
