#!/usr/bin/env python3
"""
Split gVCF reference blocks (<*>) at zero-coverage boundaries.

Reads a full-genome BedGraph (bedtools genomecov -bga) to:
  1. Identify zero-coverage regions and split/drop reference blocks accordingly.
  2. Recompute DP in the SAMPLE column for each sub-interval using the BedGraph depth.

Variant records are passed through unchanged.
"""

import sys
import gzip
import argparse
import bisect
from collections import defaultdict


# ---------------------------------------------------------------------------
# Reference FASTA loading
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# BedGraph loading
# ---------------------------------------------------------------------------

def load_bedgraph(bedgraph_file, min_coverage=1):
    """
    Load a full-genome BedGraph into memory.
    Returns:
      bgraph: dict[chrom] -> list of (start, end, depth)  (0-based half-open, sorted)
      bgraph_starts: dict[chrom] -> list of start positions (for bisect)
      low_cov_regions: dict[chrom] -> list of (start, end) where depth < min_coverage
    """
    bgraph = defaultdict(list)
    zero_regions = defaultdict(list)

    with open(bedgraph_file) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            chrom, start, end, depth = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
            bgraph[chrom].append((start, end, depth))
            if depth < min_coverage:
                zero_regions[chrom].append((start, end))

    bgraph_starts = {chrom: [e[0] for e in entries] for chrom, entries in bgraph.items()}
    return bgraph, bgraph_starts, zero_regions


# ---------------------------------------------------------------------------
# Interval arithmetic
# ---------------------------------------------------------------------------

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


def compute_depth_stats(sub_start, sub_end, bgraph_chrom, bgraph_starts_chrom):
    """
    Compute (min_dp, mean_dp) for [sub_start, sub_end] (1-based inclusive)
    using BedGraph entries (0-based half-open).
    Uses binary search to find the first overlapping entry efficiently.
    """
    q_start = sub_start - 1  # 0-based
    q_end = sub_end           # 0-based half-open
    total_bases = q_end - q_start

    # Find first entry that could overlap: last entry with start <= q_start
    idx = bisect.bisect_right(bgraph_starts_chrom, q_start) - 1
    idx = max(idx, 0)

    min_dp = None
    weighted_sum = 0

    for bg_start, bg_end, depth in bgraph_chrom[idx:]:
        if bg_end <= q_start:
            continue
        if bg_start >= q_end:
            break
        overlap = min(bg_end, q_end) - max(bg_start, q_start)
        weighted_sum += depth * overlap
        if min_dp is None or depth < min_dp:
            min_dp = depth

    mean_dp = round(weighted_sum / total_bases) if total_bases > 0 else 0
    return (min_dp if min_dp is not None else 0), mean_dp


# ---------------------------------------------------------------------------
# FORMAT/SAMPLE field updating
# ---------------------------------------------------------------------------

def update_dp(format_str, sample_str, new_dp):
    """Replace the DP value in the SAMPLE column. No-op if DP not in FORMAT."""
    fmt_fields = format_str.split(':')
    if 'DP' not in fmt_fields:
        return sample_str
    dp_idx = fmt_fields.index('DP')
    sample_fields = sample_str.split(':')
    if dp_idx < len(sample_fields):
        sample_fields[dp_idx] = str(new_dp)
    return ':'.join(sample_fields)


def update_end_in_info(info, new_end):
    fields = info.split(';')
    return ';'.join('END={}'.format(new_end) if f.startswith('END=') else f for f in fields)


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------

def process_gvcf(gvcf_file, bedgraph_file, output_file, ref_fasta=None, min_coverage=1):
    bgraph, bgraph_starts, zero_regions = load_bedgraph(bedgraph_file, min_coverage)
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
            fmt = parts[8] if len(parts) > 8 else ''
            sample = parts[9] if len(parts) > 9 else ''

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

            chrom_zeros = zero_regions.get(chrom, [])
            if not chrom_zeros:
                fout.write(line)
                continue

            intervals = covered_intervals(pos, end, chrom_zeros)

            if not intervals:
                # Entire block is zero-coverage; drop it
                continue

            if len(intervals) == 1 and intervals[0] == (pos, end):
                # No splitting needed; still update DP for accuracy
                if bgraph.get(chrom):
                    _, mean_dp = compute_depth_stats(pos, end, bgraph[chrom], bgraph_starts[chrom])
                    parts[9] = update_dp(fmt, sample, mean_dp)
                fout.write('\t'.join(parts) + '\n')
                continue

            # Emit one record per covered sub-interval with recomputed DP
            bg_chrom = bgraph.get(chrom, [])
            bg_starts_chrom = bgraph_starts.get(chrom, [])
            for sub_start, sub_end in intervals:
                new_parts = parts[:]
                new_parts[1] = str(sub_start)
                new_parts[7] = update_end_in_info(info, sub_end)
                # When the sub-interval starts at a new position, REF must reflect
                # the actual reference base there, not the base at the original block POS.
                if sub_start != pos and chrom in ref_seqs:
                    new_parts[3] = ref_seqs[chrom][sub_start - 1]  # 1-based -> 0-based
                if bg_chrom:
                    _, mean_dp = compute_depth_stats(sub_start, sub_end, bg_chrom, bg_starts_chrom)
                    new_parts[9] = update_dp(fmt, sample, mean_dp)
                fout.write('\t'.join(new_parts) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gvcf', required=True, help='Input gVCF (plain or bgzipped)')
    parser.add_argument('--bedgraph', required=True, help='Full-genome BedGraph (bedtools genomecov -bga)')
    parser.add_argument('--output', required=True, help='Output gVCF (.gz = bgzipped)')
    parser.add_argument('--ref', required=False, help='Reference FASTA; required to fix REF alleles on split sub-blocks')
    parser.add_argument('--min-coverage', required=False, type=int, default=1, dest='min_coverage',
                        help='Minimum depth to consider a region covered (default: 1)')
    args = parser.parse_args()
    process_gvcf(args.gvcf, args.bedgraph, args.output, ref_fasta=args.ref, min_coverage=args.min_coverage)


if __name__ == '__main__':
    main()
