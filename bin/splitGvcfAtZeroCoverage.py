#!/usr/bin/env python3
"""
Split gVCF reference blocks (<*>) at zero-coverage boundaries.

Reads a pre-computed zero-coverage BED file (union of all-sample zero-coverage
intervals, produced by makeMultiSampleZeroCoverageBed) to identify where
reference blocks should be split or dropped.

Variant records are passed through unchanged.
Works with both single-sample and multi-sample gVCFs.
"""

import sys
import gzip
import argparse
from collections import defaultdict


def load_fasta(fasta_file):
    """Load a FASTA file into a dict of chrom -> uppercase sequence string."""
    sequences = {}
    current_chrom = None
    chunks = []
    opener = gzip.open if fasta_file.endswith('.gz') else open
    with opener(fasta_file, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if current_chrom is not None:
                    sequences[current_chrom] = ''.join(chunks)
                current_chrom = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
    if current_chrom is not None:
        sequences[current_chrom] = ''.join(chunks)
    return sequences


def load_zero_cov_bed(bed_file):
    """
    Load a BED file of zero-coverage intervals.
    Returns: dict[chrom] -> list of (start, end)  (0-based half-open, sorted)
    """
    zero_regions = defaultdict(list)
    with open(bed_file) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            zero_regions[chrom].append((start, end))
    return {chrom: sorted(intervals) for chrom, intervals in zero_regions.items()}


def covered_intervals(pos1, end1, zero_regions_chrom):
    """
    Given a block [pos1, end1] (1-based inclusive) and zero-coverage intervals
    (0-based half-open, sorted), return list of covered (start, end) tuples
    (1-based inclusive).
    """
    block_start = pos1 - 1  # 0-based
    block_end = end1        # 0-based half-open

    covered = []
    current = block_start

    for z_start, z_end in zero_regions_chrom:
        if z_end <= block_start:
            continue
        if z_start >= block_end:
            break
        if current < z_start:
            covered.append((current + 1, z_start))  # back to 1-based
        current = max(current, z_end)

    if current < block_end:
        covered.append((current + 1, block_end))  # back to 1-based

    return covered


def update_end_in_info(info, new_end):
    fields = info.split(';')
    return ';'.join('END={}'.format(new_end) if f.startswith('END=') else f for f in fields)


def process_gvcf(gvcf_file, zero_cov_bed_file, output_file, ref_fasta=None):
    zero_regions = load_zero_cov_bed(zero_cov_bed_file)
    ref_seqs = load_fasta(ref_fasta) if ref_fasta else {}

    in_opener = gzip.open if gvcf_file.endswith('.gz') else open

    with in_opener(gvcf_file, 'rt') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue

            parts = line.rstrip('\n').split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            alt = parts[4]
            info = parts[7]

            # Pass variant records through unchanged
            if alt != '<*>':
                fout.write(line)
                continue

            # Parse END from INFO
            end = pos
            for field in info.split(';'):
                if field.startswith('END='):
                    end = int(field[4:])
                    break

            # Clamp END to chromosome length
            if chrom in ref_seqs:
                end = min(end, len(ref_seqs[chrom]))

            chrom_zeros = zero_regions.get(chrom, [])
            if not chrom_zeros:
                fout.write(line)
                continue

            intervals = covered_intervals(pos, end, chrom_zeros)

            if not intervals:
                # Entire block is zero-coverage; drop it
                continue

            if len(intervals) == 1 and intervals[0] == (pos, end):
                # No splitting needed
                fout.write(line)
                continue

            # Emit one record per covered sub-interval
            for sub_start, sub_end in intervals:
                new_parts = parts[:]
                new_parts[1] = str(sub_start)
                new_parts[7] = update_end_in_info(info, sub_end)
                if sub_start != pos and chrom in ref_seqs:
                    new_parts[3] = ref_seqs[chrom][sub_start - 1]
                fout.write('\t'.join(new_parts) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gvcf',         required=True,
                        help='Input gVCF (plain or bgzipped)')
    parser.add_argument('--zero-cov-bed', required=True, dest='zero_cov_bed',
                        help='BED file of zero-coverage intervals (union across all samples)')
    parser.add_argument('--output',       required=True,
                        help='Output gVCF')
    parser.add_argument('--ref',          required=False,
                        help='Reference FASTA; fixes REF alleles on split sub-blocks')
    args = parser.parse_args()
    process_gvcf(args.gvcf, args.zero_cov_bed, args.output, ref_fasta=args.ref)


if __name__ == '__main__':
    main()
