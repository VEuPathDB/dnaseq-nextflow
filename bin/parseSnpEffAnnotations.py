#!/usr/bin/env python3
"""Parse SnpEff ANN INFO field from annotated VCF and write snpeff.dat.

Output columns: location, seq_id, allele, transcript_id, snpeff_impact, snpeff_effect
One row per unique (location, seq_id, allele, transcript_id) combination.
"""
import argparse
import gzip


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf",    required=True, help="SnpEff-annotated VCF (.vcf or .vcf.gz)")
    p.add_argument("--output", required=True, help="Output path for snpeff.dat")
    return p.parse_args()


def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    args = parse_args()
    seen = set()

    with open_vcf(args.vcf) as vcf_fh, open(args.output, "w") as out:
        out.write("location\tseq_id\tallele\ttranscript_id\tsnpeff_impact\tsnpeff_effect\n")
        for line in vcf_fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            seq_id   = fields[0]
            location = fields[1]
            info     = fields[7]

            ann_value = None
            for token in info.split(";"):
                if token.startswith("ANN="):
                    ann_value = token[4:]
                    break
            if ann_value is None:
                continue

            for entry in ann_value.split(","):
                parts = [p.strip() for p in entry.split("|")]
                if len(parts) < 7:
                    continue
                allele        = parts[0]
                effect        = parts[1]
                impact        = parts[2]
                transcript_id = parts[6]
                if not allele or not transcript_id or not impact or not effect:
                    continue
                # One row per (location, seq_id, allele, transcript_id). When SnpEff emits
                # multiple effects for the same transcript+allele, first one wins.
                # SnpEff orders ANN entries by impact severity (HIGH first), so first-wins
                # preserves the most severe annotation for each transcript.
                key = (location, seq_id, allele, transcript_id)
                if key in seen:
                    continue
                seen.add(key)
                out.write(f"{location}\t{seq_id}\t{allele}\t{transcript_id}\t{impact}\t{effect}\n")


if __name__ == "__main__":
    main()
