#!/usr/bin/env python3
"""Convert a FreeBayes VCF to VarScan-compatible format.

Makes two changes so that tools expecting VarScan output
(e.g. makeHeterozygosityPlot.py) can parse the file:

1. Adds RD FORMAT field (copy of FreeBayes RO = reference allele depth).
2. Rewrites the AD field from FreeBayes comma-separated allele depths
   (ref,alt1,alt2,...) to a single integer (sum of alt depths), matching
   VarScan's AD which is alternate-only depth.
"""

import argparse
import gzip
import re
import sys


def open_vcf(path, mode='rt'):
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vcfFile', required=True, help='Input FreeBayes VCF (.vcf or .vcf.gz)')
    args = parser.parse_args()

    fout = sys.stdout

    with open_vcf(args.vcfFile) as fin:
        for line in fin:
            # Rewrite AD header: Number=R means "one value per allele" (a list).
            # Change to Number=1 so parsers treat it as a scalar alt depth.
            if line.startswith('##FORMAT=<ID=AD,'):
                line = re.sub(r'Number=[^,]+', 'Number=1', line, count=1)
                fout.write(line)
                continue

            # Duplicate the RO FORMAT header line as RD
            if line.startswith('##FORMAT=<ID=RO,'):
                fout.write(line)
                fout.write(line.replace('ID=RO,', 'ID=RD,', 1))
                continue

            if line.startswith('#'):
                fout.write(line)
                continue

            fields = line.rstrip('\n').split('\t')
            fmt = fields[8].split(':')

            ro_idx = fmt.index('RO') if 'RO' in fmt else None
            ad_idx = fmt.index('AD') if 'AD' in fmt else None

            if ro_idx is not None:
                fields[8] += ':RD'

            for i in range(9, len(fields)):
                sample = fields[i].split(':')

                # Rewrite AD: FreeBayes uses ref,alt1,alt2,... — collapse alts to
                # a single integer to match VarScan's single-value AD field.
                if ad_idx is not None and ad_idx < len(sample):
                    ad_parts = sample[ad_idx].split(',')
                    if len(ad_parts) > 1:
                        alt_depths = [int(x) for x in ad_parts[1:] if x != '.']
                        sample[ad_idx] = str(sum(alt_depths)) if alt_depths else '0'

                # Append RD = RO value
                if ro_idx is not None:
                    ro_val = sample[ro_idx] if ro_idx < len(sample) else '.'
                    sample.append(ro_val)

                fields[i] = ':'.join(sample)

            fout.write('\t'.join(fields) + '\n')

        fout.flush()


if __name__ == '__main__':
    main()
