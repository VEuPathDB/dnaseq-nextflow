#!/usr/bin/env python3
"""
Reads VCF from stdin, writes to stdout.

After bcftools norm -a decomposes a multi-allelic site into per-allele rows,
each row carries * as a placeholder for the alleles contributed by the other rows.
This script groups rows sharing the same (CHROM, POS, REF) and recombines them
into a single multi-allelic record with a corrected GT.

For a single row with * in ALT:
  - Remaps GT allele indices that point to * back to the non-* allele
  - Drops records where * is the only ALT
"""

import sys


def _gt_sep(gt_str):
    return '|' if '|' in gt_str else '/'


def _fix_single(fields):
    """Handle a single row with * in ALT. Returns updated fields list or None to drop."""
    alts = fields[4].split(',')

    if '*' not in alts:
        return fields

    if alts == ['*']:
        return None

    star_alt_pos = alts.index('*')
    star_gt_idx  = star_alt_pos + 1  # 1-based GT index for *

    new_alts = []
    old_to_new = {0: 0}
    new_idx = 1
    for i, a in enumerate(alts):
        if a == '*':
            old_to_new[i + 1] = 1
        else:
            old_to_new[i + 1] = new_idx
            new_alts.append(a)
            new_idx += 1

    fields[4] = ','.join(new_alts)

    n_alts  = len(alts)       # original ALT count (includes *)
    n_total = n_alts + 1      # total alleles including REF

    # For Number=G (diploid): precompute which genotype indices to keep
    # (those where neither allele is *)
    gl_keep = [
        j * (j + 1) // 2 + i
        for j in range(n_total)
        for i in range(j + 1)
        if i != star_gt_idx and j != star_gt_idx
    ]

    fmt    = fields[8].split(':')
    gt_pos = fmt.index('GT') if 'GT' in fmt else 0

    for col in range(9, len(fields)):
        parts     = fields[col].split(':')
        new_parts = []
        for fi, tag in enumerate(fmt):
            val  = parts[fi] if fi < len(parts) else '.'
            if tag == 'GT':
                sep = _gt_sep(val)
                new_parts.append(sep.join(
                    '.' if a == '.' else str(old_to_new.get(int(a), int(a)))
                    for a in val.split(sep)
                ))
            else:
                vals = val.split(',')
                if len(vals) == n_alts + 1:        # Number=R (REF + each ALT)
                    vals.pop(star_alt_pos + 1)
                elif len(vals) == n_alts:           # Number=A (one per ALT)
                    vals.pop(star_alt_pos)
                elif len(vals) == n_total * (n_total + 1) // 2:  # Number=G diploid
                    vals = [vals[k] for k in gl_keep if k < len(vals)]
                new_parts.append(','.join(vals))
        fields[col] = ':'.join(new_parts)

    return fields


def _recombine(group_fields):
    """
    Recombine multiple rows at the same (CHROM, POS, REF) into one row.

    Row i contributes one real (non-*) ALT allele, which becomes combined ALT
    index i+1.  The combined GT is built by taking one allele index per row —
    the one pointing to that row's real allele — and mapping it to i+1.

    For per-allele FORMAT fields (Number=R, Number=A), * placeholder values are
    ignored and real values are pulled from each row's own real allele position.
    bcftools -a guarantees at most 2 rows per site, each with exactly one real ALT.
    """
    real_alleles = []
    real_alt_indices = []  # 0-based index in each row's ALT list for the non-* allele
    for fields in group_fields:
        alts = fields[4].split(',')
        real = [i for i, a in enumerate(alts) if a != '*']
        real_alleles.extend(alts[i] for i in real)
        real_alt_indices.append(real[0] if real else None)

    if not real_alleles:
        return None

    template = list(group_fields[0])
    template[4] = ','.join(real_alleles)

    fmt = template[8].split(':')
    gt_pos = fmt.index('GT') if 'GT' in fmt else 0

    first_gt = group_fields[0][9].split(':')[gt_pos]
    sep = _gt_sep(first_gt)

    combined_gt_parts = []
    for i, fields in enumerate(group_fields):
        alts = fields[4].split(',')
        fmt_i = fields[8].split(':')
        gp_i = fmt_i.index('GT') if 'GT' in fmt_i else 0
        gt = fields[9].split(':')[gp_i]
        before = len(combined_gt_parts)
        for allele in gt.split(_gt_sep(gt)):
            if allele == '.':
                continue
            idx = int(allele)
            if idx == 0:
                combined_gt_parts.append('0')
                break
            if idx <= len(alts) and alts[idx - 1] != '*':
                combined_gt_parts.append(str(i + 1))
                break
        # All GT alleles pointed to * so nothing was appended. This is an artifact of
        # bcftools norm -a decomposing a complex variant (e.g. MNP + indel at the same
        # position) where neither decomposed row carries the real allele in its GT.
        # Fall back to assigning this row's allele positionally (1/2 for a 2-row group).
        if len(combined_gt_parts) == before:
            combined_gt_parts.append(str(i + 1))

    n_alts_row0 = len(group_fields[0][4].split(','))

    for col in range(9, len(template)):
        row0_parts = group_fields[0][col].split(':')
        new_parts = []
        for fi, tag in enumerate(fmt):
            if tag == 'GT':
                new_parts.append(sep.join(combined_gt_parts))
                continue

            val0 = (row0_parts[fi] if fi < len(row0_parts) else '.').split(',')

            if len(val0) == n_alts_row0 + 1:  # Number=R: REF + one per ALT
                ref_val = val0[0]
                alt_vals = []
                for ri, rfields in enumerate(group_fields):
                    n_alts_ri = len(rfields[4].split(','))
                    rparts = rfields[col].split(':')
                    rv = (rparts[fi] if fi < len(rparts) else '.').split(',')
                    ridx = real_alt_indices[ri]
                    alt_vals.append(rv[ridx + 1] if ridx is not None and len(rv) == n_alts_ri + 1 else '.')
                new_parts.append(','.join([ref_val] + alt_vals))

            elif len(val0) == n_alts_row0:  # Number=A: one per ALT
                alt_vals = []
                for ri, rfields in enumerate(group_fields):
                    n_alts_ri = len(rfields[4].split(','))
                    rparts = rfields[col].split(':')
                    rv = (rparts[fi] if fi < len(rparts) else '.').split(',')
                    ridx = real_alt_indices[ri]
                    alt_vals.append(rv[ridx] if ridx is not None and len(rv) == n_alts_ri else '.')
                new_parts.append(','.join(alt_vals))

            else:  # Number=1, Number=G, etc. — keep row 0's value
                new_parts.append(row0_parts[fi] if fi < len(row0_parts) else '.')

        template[col] = ':'.join(new_parts)

    return template


def _flush(group):
    """Process a buffered group of rows and return output lines."""
    if not group:
        return []

    parsed = [row.rstrip('\n').split('\t') for row in group]
    all_have_star = all('*' in f[4].split(',') for f in parsed if len(f) >= 5)

    if len(parsed) > 1 and all_have_star:
        result = _recombine(parsed)
        return ['\t'.join(result) + '\n'] if result else []

    out = []
    for fields in parsed:
        fixed = _fix_single(fields)
        if fixed:
            out.append('\t'.join(fixed) + '\n')
    return out


def main():
    buffer = []
    current_key = None

    for line in sys.stdin:
        if line.startswith('#'):
            for l in _flush(buffer):
                sys.stdout.write(l)
            buffer = []
            current_key = None
            sys.stdout.write(line)
            continue

        fields = line.rstrip('\n').split('\t')
        if len(fields) < 5:
            sys.stdout.write(line)
            continue

        key = (fields[0], fields[1], fields[3])  # CHROM, POS, REF

        if key != current_key:
            for l in _flush(buffer):
                sys.stdout.write(l)
            buffer = [line]
            current_key = key
        else:
            buffer.append(line)

    for l in _flush(buffer):
        sys.stdout.write(l)


if __name__ == '__main__':
    main()
