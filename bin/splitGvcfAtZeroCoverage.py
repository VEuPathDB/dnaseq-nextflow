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
import os
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


def load_per_sample_beds(bed_paths):
    """
    Load per-sample zero-coverage BEDs named <samplename>.persample.bed.
    Returns dict[sample_name] -> dict[chrom] -> [(start, end)]
    """
    per_sample = {}
    for bed_path in (bed_paths or []):
        basename = os.path.basename(bed_path)
        if not basename.endswith('.persample.bed'):
            continue
        sample_name = basename[:-len('.persample.bed')]
        per_sample[sample_name] = load_zero_cov_bed(bed_path)
    return per_sample


def parse_sample_names(gvcf_file):
    """Extract ordered sample names from gVCF #CHROM header line."""
    opener = gzip.open if gvcf_file.endswith('.gz') else open
    with opener(gvcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                cols = line.rstrip('\n').split('\t')
                return cols[9:]
    return []


def all_sub_intervals(pos1, end1, zero_regions_chrom):
    """
    Split [pos1, end1] (1-based inclusive) at zero-coverage boundaries.
    Returns ALL sub-intervals in order — both covered and zero-coverage — as
    (start, end) tuples (1-based inclusive).
    """
    block_start = pos1 - 1  # 0-based
    block_end = end1        # 0-based half-open

    boundaries = {block_start, block_end}
    for z_start, z_end in zero_regions_chrom:
        if z_end <= block_start or z_start >= block_end:
            continue
        boundaries.add(max(z_start, block_start))
        boundaries.add(min(z_end, block_end))

    sorted_bounds = sorted(boundaries)
    return [(sorted_bounds[i] + 1, sorted_bounds[i + 1])
            for i in range(len(sorted_bounds) - 1)]


def is_all_zero(sub_start, sub_end, all_zero_chrom):
    """
    Return True if [sub_start, sub_end] (1-based inclusive) is entirely covered
    by at least one interval in all_zero_chrom (0-based half-open, sorted).
    """
    b_start = sub_start - 1  # 0-based
    b_end = sub_end          # 0-based half-open
    for z_start, z_end in all_zero_chrom:
        if z_start <= b_start and z_end >= b_end:
            return True
        if z_start >= b_end:
            break
    return False


def update_end_in_info(info, new_end):
    fields = info.split(';')
    return ';'.join('END={}'.format(new_end) if f.startswith('END=') else f for f in fields)


def zero_sample_dp(format_str, sample_cols, sample_names, per_sample_chrom, sub_start, sub_end):
    """
    Return updated sample columns with DP (and MIN_DP) set to 0 for any sample
    whose per-sample zero-coverage BED covers [sub_start, sub_end] entirely.
    Leaves '.' values and already-zero values untouched.
    """
    fmt_keys = format_str.split(':')
    dp_idx     = fmt_keys.index('DP')     if 'DP'     in fmt_keys else None
    min_dp_idx = fmt_keys.index('MIN_DP') if 'MIN_DP' in fmt_keys else None
    if dp_idx is None:
        return sample_cols

    new_cols = list(sample_cols)
    for i, sample_name in enumerate(sample_names):
        chrom_regions = per_sample_chrom.get(sample_name, [])
        if not chrom_regions:
            continue
        if not is_all_zero(sub_start, sub_end, chrom_regions):
            continue
        vals = new_cols[i].split(':')
        if dp_idx < len(vals) and vals[dp_idx] not in ('.', '0'):
            vals[dp_idx] = '0'
        if min_dp_idx is not None and min_dp_idx < len(vals) and vals[min_dp_idx] not in ('.', '0'):
            vals[min_dp_idx] = '0'
        new_cols[i] = ':'.join(vals)
    return new_cols


def process_gvcf(gvcf_file, union_zero_bed_file, all_zero_bed_file, output_file,
                 ref_fasta=None, per_sample_beds=None):
    union_regions    = load_zero_cov_bed(union_zero_bed_file)
    all_zero_regions = load_zero_cov_bed(all_zero_bed_file)
    ref_seqs = load_fasta(ref_fasta) if ref_fasta else {}

    per_sample_zero = load_per_sample_beds(per_sample_beds)
    sample_names    = parse_sample_names(gvcf_file)

    in_opener = gzip.open if gvcf_file.endswith('.gz') else open

    with in_opener(gvcf_file, 'rt') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue

            parts = line.rstrip('\n').split('\t')
            chrom      = parts[0]
            pos        = int(parts[1])
            alt        = parts[4]
            info       = parts[7]
            format_str = parts[8]
            sample_cols = parts[9:]

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

            union_chrom    = union_regions.get(chrom, [])
            all_zero_chrom = all_zero_regions.get(chrom, [])
            per_sample_chrom = {s: per_sample_zero[s].get(chrom, []) for s in per_sample_zero}

            if not union_chrom:
                # No split boundaries at all; only drop if entirely all-zero
                if not is_all_zero(pos, end, all_zero_chrom):
                    fout.write(line)
                continue

            # Split at union boundaries, then drop only all-zero sub-intervals
            sub_intervals = all_sub_intervals(pos, end, union_chrom)

            for sub_start, sub_end in sub_intervals:
                if is_all_zero(sub_start, sub_end, all_zero_chrom):
                    continue
                new_parts = parts[:]
                new_parts[1] = str(sub_start)
                new_parts[7] = update_end_in_info(info, sub_end)
                if sub_start != pos and chrom in ref_seqs:
                    new_parts[3] = ref_seqs[chrom][sub_start - 1]
                updated_samples = zero_sample_dp(
                    format_str, sample_cols, sample_names, per_sample_chrom, sub_start, sub_end
                )
                new_parts[9:] = updated_samples
                fout.write('\t'.join(new_parts) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gvcf',             required=True,
                        help='Input gVCF (plain or bgzipped)')
    parser.add_argument('--union-zero-bed',   required=True, dest='union_zero_bed',
                        help='BED of positions zero in ANY present sample (split boundaries)')
    parser.add_argument('--all-zero-bed',     required=True, dest='all_zero_bed',
                        help='BED of positions zero in ALL present samples (drop decisions)')
    parser.add_argument('--per-sample-beds',  required=False, nargs='*', dest='per_sample_beds',
                        help='Per-sample zero-coverage BEDs named <samplename>.persample.bed')
    parser.add_argument('--output',           required=True,
                        help='Output gVCF')
    parser.add_argument('--ref',              required=False,
                        help='Reference FASTA; fixes REF alleles on split sub-blocks')
    args = parser.parse_args()
    process_gvcf(args.gvcf, args.union_zero_bed, args.all_zero_bed, args.output,
                 ref_fasta=args.ref, per_sample_beds=args.per_sample_beds)


if __name__ == '__main__':
    main()
