#!/usr/bin/env python3
"""
Reads VCF from stdin, writes to stdout.
For records with * (spanning deletion) in the ALT field:
  - Remaps GT allele indices that point to * back to the non-* allele
  - Removes * from the ALT field
  - Drops records where * is the only ALT
"""

import sys


def fix_record(line):
    if line.startswith('#'):
        return line

    fields = line.rstrip('\n').split('\t')
    if len(fields) < 10:
        return line

    alts = fields[4].split(',')

    if '*' not in alts:
        return line

    # Drop records where * is the only ALT
    if alts == ['*']:
        return None

    # Index of * in the ALT list (0-based), and its GT allele index (1-based, REF=0)
    star_alt_pos = alts.index('*')
    star_gt_idx  = star_alt_pos + 1

    # Build new ALT without *, and a mapping from old GT index -> new GT index
    new_alts = []
    old_to_new = {0: 0}  # REF stays 0
    new_idx = 1
    for i, a in enumerate(alts):
        old_gt_idx = i + 1
        if a == '*':
            old_to_new[old_gt_idx] = 1  # remap * to non-* allele (always index 1 post-removal)
        else:
            old_to_new[old_gt_idx] = new_idx
            new_alts.append(a)
            new_idx += 1

    fields[4] = ','.join(new_alts)

    # Fix GT in each sample column
    fmt = fields[8].split(':')
    gt_pos = fmt.index('GT') if 'GT' in fmt else 0

    for col in range(9, len(fields)):
        parts = fields[col].split(':')
        gt = parts[gt_pos]
        sep = '|' if '|' in gt else '/'
        alleles = gt.split(sep)
        new_alleles = []
        for a in alleles:
            if a == '.':
                new_alleles.append('.')
            else:
                new_alleles.append(str(old_to_new.get(int(a), int(a))))
        parts[gt_pos] = sep.join(new_alleles)
        fields[col] = ':'.join(parts)

    return '\t'.join(fields) + '\n'


def main():
    for line in sys.stdin:
        result = fix_record(line)
        if result is not None:
            sys.stdout.write(result)


if __name__ == '__main__':
    main()
