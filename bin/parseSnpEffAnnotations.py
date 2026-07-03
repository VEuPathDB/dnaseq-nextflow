#!/usr/bin/env python3
"""Parse coding annotations from an annotated VCF and write snpeff.dat.

Two annotation sources are combined into one file, distinguished by the
`source` column:
  - snpeff:       SnpEff's ANN INFO field (functional annotation)
  - product_call: this pipeline's hand-rolled coding calls carried in the
                  CANN INFO field (see processSequenceVariations.jl)

Output columns:
  location, seq_id, allele, transcript_id, impact, effect, source
One row per unique (location, seq_id, allele, transcript_id, effect, source).

Compound effects (SnpEff joins co-occurring consequences with '&', e.g.
"missense_variant&splice_region_variant") are split into one row per effect so
each row carries a single effect and a single impact tier — clean for equality
filters in the consuming webapp. SnpEff reports one impact per entry, so its
split rows replicate that entry-level impact; product_call rows get the exact
tier for each component.

The consumer loads this into a relational DB; multiple rows per position are
expected and queries match on location.
"""
import argparse
import gzip


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf",    required=True, help="Annotated VCF (.vcf or .vcf.gz)")
    p.add_argument("--output", required=True, help="Output path for snpeff.dat")
    return p.parse_args()


def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


# ---------------------------------------------------------------------------
# CANN (product_call) effect vocabulary
#
# Maps the hand-rolled effect terms emitted by processSequenceVariations.jl's
# build_cann_string onto SnpEff/Sequence-Ontology effect terms and an impact
# tier, so both sources speak the same vocabulary in the shared file.
# ---------------------------------------------------------------------------

CANN_EFFECT_MAP = {
    "synonymous":            ("synonymous_variant",    "LOW"),
    "missense":              ("missense_variant",      "MODERATE"),
    "nonsense":              ("stop_gained",           "HIGH"),
    "frameshift":            ("frameshift_variant",    "HIGH"),
    # SnpEff always classifies these as conservative_/disruptive_ after
    # renormalizing indels to a canonical position within repeat regions.
    # Julia doesn't do that renormalization, so its codon-boundary math can't
    # reliably tell conservative from disruptive apart (confirmed: it disagrees
    # with SnpEff on most repeat-region indels in test data). Emit a distinct
    # term rather than mislabel severity or falsely imply parity with SnpEff's
    # classified calls.
    "inframe_insertion":     ("inframe_insertion_unnormalized", "MODERATE"),
    "inframe_deletion":      ("inframe_deletion_unnormalized",  "MODERATE"),
    # Positional consequence flag, not a direct sequence change.
    "downstream_frameshift": ("downstream_frameshift", "MODIFIER"),
}

# Effect ".": Julia could not resolve the codon. The codon field disambiguates:
#   - contains X  -> a het indel masked the consensus (X-masking rule); queryable
#                    as its own term.
#   - contains N  -> missing strain sequence; consequence genuinely unknown.
CANN_HET_INDEL  = ("het_indel", "MODIFIER")
CANN_UNRESOLVED = ("coding_sequence_variant", "MODIFIER")


def map_cann_effect(effect, codon):
    """Map a Julia effect string to a list of (so_effect, impact_tier) pairs,
    one per '&'-joined component. An unresolved effect (".") is classified from
    the codon: an X means a het indel, otherwise consequence is unknown."""
    if not effect or effect == ".":
        return [CANN_HET_INDEL if "X" in codon.upper() else CANN_UNRESOLVED]
    return [CANN_EFFECT_MAP.get(comp, CANN_UNRESOLVED) for comp in effect.split("&")]


def parse_ann_rows(info):
    """Yield (allele, transcript_id, impact, effect) tuples from the ANN field.

    Compound effects are split on '&'; each component row inherits the entry's
    single SnpEff impact.
    """
    ann_value = None
    for token in info.split(";"):
        if token.startswith("ANN="):
            ann_value = token[4:]
            break
    if ann_value is None:
        return

    for entry in ann_value.split(","):
        parts = [p.strip() for p in entry.split("|")]
        if len(parts) < 7:
            continue
        allele        = parts[0]
        effect        = parts[1]
        impact        = parts[2]
        feature_type  = parts[5]
        # parts[6] is Feature_ID — a real transcript ID only when Feature_Type is
        # "transcript"; for intergenic_region SnpEff puts a gene-boundary name there.
        transcript_id = parts[6] if feature_type == "transcript" else ""
        if not allele or not impact or not effect:
            continue
        if feature_type == "transcript" and not transcript_id:
            continue
        for component in effect.split("&"):
            yield (allele, transcript_id, impact, component)


def parse_cann_rows(alt, info):
    """Yield (allele, transcript_id, impact, effect) tuples from the CANN field.

    Each output VCF line carries a single ALT, so every k-prefixed CANN entry
    on the line pertains to that ALT. r-prefixed entries describe the reference
    allele and are skipped. Compound effects are split, each with its own tier.
    CANN entry format: key|codon|aa|effect|transcript_id|pos_in_cds|pos_in_codon
    """
    cann_value = None
    for token in info.split(";"):
        if token.startswith("CANN="):
            cann_value = token[5:]
            break
    if cann_value is None or cann_value == ".":
        return

    for entry in cann_value.split(","):
        parts = entry.split("|")
        if len(parts) < 5:
            continue
        key           = parts[0]
        codon         = parts[1]
        effect        = parts[3]
        transcript_id = parts[4]
        if not key.startswith("k"):   # skip reference (r-prefixed) entries
            continue
        if not transcript_id:
            continue
        for so_effect, impact in map_cann_effect(effect, codon):
            yield (alt, transcript_id, impact, so_effect)


def main():
    args = parse_args()
    seen = set()

    with open_vcf(args.vcf) as vcf_fh, open(args.output, "w") as out:
        out.write("location\tseq_id\tallele\ttranscript_id\timpact\teffect\tsource\n")
        for line in vcf_fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            seq_id   = fields[0]
            location = fields[1]
            alt      = fields[4]
            info     = fields[7]

            sourced_rows = [
                ("snpeff",       row) for row in parse_ann_rows(info)
            ] + [
                ("product_call", row) for row in parse_cann_rows(alt, info)
            ]

            for source, (allele, transcript_id, impact, effect) in sourced_rows:
                # All distinct effects for a transcript+allele are kept — `effect`
                # is part of the key. Only an exact-duplicate row (same effect and
                # source) is collapsed, since it would be byte-identical output.
                key = (location, seq_id, allele, transcript_id, effect, source)
                if key in seen:
                    continue
                seen.add(key)
                out.write(f"{location}\t{seq_id}\t{allele}\t{transcript_id}\t{impact}\t{effect}\t{source}\n")


if __name__ == "__main__":
    main()
