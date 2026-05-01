#!/usr/bin/env python3
import argparse
import os
import subprocess

try:
    from cyvcf2 import VCF
except ImportError:
    VCF = None

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


def build_consensus(chrom_name, chrom_len, ref_seq, vcf, min_coverage, sample_idx):
    """
    Walk through g.vcf records for one chromosome and assemble the consensus
    for the sample at sample_idx.
    """
    segments = []
    ref_pos = 0

    for v in vcf(chrom_name):
        pos = v.POS - 1  # VCF 1-based → 0-based

        if pos > ref_pos:
            segments.append('N' * (pos - ref_pos))
            ref_pos = pos

        dp_arr = v.format('DP')
        if dp_arr is None:
            dp = 0
        else:
            raw = int(dp_arr[sample_idx][0])
            dp = 0 if raw < 0 else raw  # cyvcf2 uses negative sentinel for missing ('.')

        # REF block
        if not v.ALT or all(a == '.' or a.startswith('<') for a in v.ALT):
            end = (v.INFO.get('END') or v.POS) - 1  # 0-based inclusive
            span = end - pos + 1
            if dp >= min_coverage:
                segments.append(ref_seq[pos:end + 1])
            else:
                segments.append('N' * span)
            ref_pos = end + 1
            continue

        # Variant record
        if dp < min_coverage:
            segments.append('N' * len(v.REF))
            ref_pos = pos + len(v.REF)
            continue

        if 'GT' not in v.FORMAT:
            segments.append('N' * len(v.REF))
            ref_pos = pos + len(v.REF)
            continue

        gt_str = v.gt_bases[sample_idx]
        if '.' in gt_str:
            segments.append('N' * len(v.REF))
            ref_pos = pos + len(v.REF)
            continue

        alleles = list(dict.fromkeys(gt_str.replace('|', '/').split('/')))

        if all(len(a) == 1 for a in alleles):
            base = IUPAC.get(frozenset(alleles), 'N')
            segments.append(base)
            ref_pos = pos + 1

        elif len(alleles) == 1:
            segments.append(alleles[0])
            ref_pos = pos + len(v.REF)

        else:
            if v.REF in alleles:
                segments.append(v.REF)
            else:
                segments.append('X' * len(v.REF))
            ref_pos = pos + len(v.REF)

    if ref_pos < chrom_len:
        segments.append('N' * (chrom_len - ref_pos))

    return ''.join(segments)


def write_fasta(out_fh, name, seq, line_len=60):
    out_fh.write(f'>{name}\n')
    for i in range(0, len(seq), line_len):
        out_fh.write(seq[i:i + line_len] + '\n')


def main():
    if VCF is None:
        raise RuntimeError('cyvcf2 is required to run this script.')

    parser = argparse.ArgumentParser(
        description='Build per-sample consensus FASTAs from a multi-sample g.vcf.')
    parser.add_argument('-g', '--gvcf',          required=True,
                        help='g.vcf.gz with tabix index (.tbi)')
    parser.add_argument('-r', '--ref',           required=True,
                        help='Reference FASTA (indexed with samtools faidx)')
    parser.add_argument('-f', '--fai',           required=True,
                        help='Reference FASTA .fai index')
    parser.add_argument('-mc', '--min-coverage', required=True, type=int,
                        dest='min_coverage')
    parser.add_argument('-o', '--output-dir',    required=False, default='.',
                        dest='output_dir',
                        help='Directory for output FASTA files (default: current directory)')
    args = parser.parse_args()

    chroms = []
    with open(args.fai) as fh:
        for line in fh:
            parts = line.split('\t')
            chroms.append((parts[0], int(parts[1])))

    # Read sample names from header, then iterate once per sample
    with VCF(args.gvcf) as hdr:
        sample_names = list(hdr.samples)

    os.makedirs(args.output_dir, exist_ok=True)

    for sample_idx, sample_name in enumerate(sample_names):
        vcf = VCF(args.gvcf)
        out_path = os.path.join(args.output_dir, f'{sample_name}_consensus.fa')
        with open(out_path, 'w') as out:
            for chrom_name, chrom_len in chroms:
                ref_seq = get_chrom_seq(args.ref, chrom_name)
                seq = build_consensus(chrom_name, chrom_len, ref_seq, vcf,
                                      args.min_coverage, sample_idx)
                write_fasta(out, chrom_name, seq)


if __name__ == '__main__':
    main()
