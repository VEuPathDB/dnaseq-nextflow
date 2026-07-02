"""Tests for bin/parseSnpEffAnnotations.py"""
import gzip
import os
import subprocess
import tempfile

import pytest

SCRIPT = os.path.join(os.path.dirname(__file__), "../../bin/parseSnpEffAnnotations.py")


def run_script(vcf_content: str, compressed: bool = True, include_header: bool = False) -> list[str]:
    """Write VCF content to a temp file, run script, return output lines.

    Header is dropped unless include_header=True.
    """
    with tempfile.TemporaryDirectory() as tmp:
        vcf_path = os.path.join(tmp, "test.vcf.gz" if compressed else "test.vcf")
        out_path  = os.path.join(tmp, "snpeff.dat")
        if compressed:
            with gzip.open(vcf_path, 'wt') as f:
                f.write(vcf_content)
        else:
            with open(vcf_path, 'w') as f:
                f.write(vcf_content)
        subprocess.run(
            ["python3", SCRIPT, "--vcf", vcf_path, "--output", out_path],
            check=True
        )
        with open(out_path) as f:
            lines = f.read().strip().split('\n')
        return lines if include_header else lines[1:]


MINIMAL_HEADER = "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

ANN_MODERATE    = "T|missense_variant|MODERATE|geneA|geneA_id|transcript|tx1|protein_coding|1/3|c.100A>T|p.Ser34Cys|100/900|100/900|34/300|.|"
ANN_HIGH        = "A|stop_gained|HIGH|geneA|geneA_id|transcript|tx1|protein_coding|1/3|c.50G>A|p.Trp17*|50/900|50/900|17/300|.|"
ANN_INTERGENIC  = "G|intergenic_region|MODIFIER|CHR_START-LmjF.01.0010|CHR_START-LmjF.01.0010|intergenic_region|CHR_START-LmjF.01.0010|||n.233C>G||||||"


def test_parses_single_ann_entry():
    vcf = MINIMAL_HEADER + f"LmjF.01\t3745\t.\tC\tT\t.\t.\tANN={ANN_MODERATE}\n"
    rows = run_script(vcf)
    assert len(rows) == 1
    fields = rows[0].split('\t')
    assert fields[0] == "3745"           # location
    assert fields[1] == "LmjF.01"        # seq_id
    assert fields[2] == "T"              # allele
    assert fields[3] == "tx1"            # transcript_id
    assert fields[4] == "MODERATE"       # snpeff_impact
    assert fields[5] == "missense_variant"  # snpeff_effect


def test_parses_multiple_ann_entries_same_position():
    vcf = MINIMAL_HEADER + f"LmjF.01\t3745\t.\tG\tT,A\t.\t.\tANN={ANN_MODERATE},{ANN_HIGH}\n"
    rows = run_script(vcf)
    assert len(rows) == 2
    alleles = {r.split('\t')[2] for r in rows}
    impacts = {r.split('\t')[4] for r in rows}
    assert alleles == {"T", "A"}
    assert impacts == {"MODERATE", "HIGH"}


def test_skips_entries_with_empty_transcript():
    no_tx = "T|missense_variant|MODERATE|geneA|geneA_id|transcript||protein_coding|1/3|.|.|.|.|.|.|"
    vcf   = MINIMAL_HEADER + f"LmjF.01\t100\t.\tC\tT\t.\t.\tANN={no_tx}\n"
    rows  = [r for r in run_script(vcf) if r]
    assert rows == []   # no data rows (empty transcript is skipped)


def test_deduplicates_same_allele_transcript_pair():
    dup = f"{ANN_MODERATE},{ANN_MODERATE}"
    vcf = MINIMAL_HEADER + f"LmjF.01\t200\t.\tC\tT\t.\t.\tANN={dup}\n"
    rows = [r for r in run_script(vcf) if r]
    assert len(rows) == 1


def test_handles_uncompressed_vcf():
    vcf  = MINIMAL_HEADER + f"LmjF.01\t500\t.\tC\tT\t.\t.\tANN={ANN_MODERATE}\n"
    rows = run_script(vcf, compressed=False)
    assert len(rows) == 1
    assert rows[0].split('\t')[4] == "MODERATE"


def test_intergenic_written_with_empty_transcript_id():
    vcf  = MINIMAL_HEADER + f"LmjF.01\t233\t.\tC\tG\t.\t.\tANN={ANN_INTERGENIC}\n"
    rows = [r for r in run_script(vcf) if r]
    assert len(rows) == 1
    fields = rows[0].split('\t')
    assert fields[0] == "233"
    assert fields[3] == ""                  # transcript_id must be empty
    assert fields[4] == "MODIFIER"
    assert fields[5] == "intergenic_region"


def test_skips_records_without_ann():
    vcf = MINIMAL_HEADER + "LmjF.01\t999\t.\tC\tT\t.\t.\tDP=30\n"
    rows = [r for r in run_script(vcf) if r]
    assert rows == []


# ---------------------------------------------------------------------------
# source column + hand-rolled product calls from the CANN INFO field
# CANN entry format: key|codon|aa|effect|transcript_id|pos_in_cds|pos_in_codon
# ---------------------------------------------------------------------------

def test_ann_row_has_source_snpeff():
    vcf = MINIMAL_HEADER + f"LmjF.01\t3745\t.\tC\tT\t.\t.\tANN={ANN_MODERATE}\n"
    rows = run_script(vcf)
    assert len(rows) == 1
    fields = rows[0].split('\t')
    assert fields[6] == "snpeff"      # new source column


def test_cann_missense_becomes_product_call_row():
    vcf = MINIMAL_HEADER + "LmjF.01\t100\t.\tC\tT\t.\t.\tCANN=k0|ATG|M|missense|tx1|100|1\n"
    rows = [r for r in run_script(vcf) if r]
    assert len(rows) == 1
    f = rows[0].split('\t')
    assert f[0] == "100"                 # location
    assert f[1] == "LmjF.01"             # seq_id
    assert f[2] == "T"                   # allele (VCF ALT)
    assert f[3] == "tx1"                 # transcript_id
    assert f[4] == "MODERATE"            # impact tier
    assert f[5] == "missense_variant"    # SnpEff-aligned effect term
    assert f[6] == "product_call"        # source


def test_cann_effect_and_impact_mapping():
    cases = {
        "synonymous":           ("synonymous_variant", "LOW"),
        "missense":             ("missense_variant",   "MODERATE"),
        "nonsense":             ("stop_gained",        "HIGH"),
        "frameshift":           ("frameshift_variant", "HIGH"),
        "inframe_insertion":    ("inframe_insertion",  "MODERATE"),
        "inframe_deletion":     ("inframe_deletion",   "MODERATE"),
        "downstream_frameshift":("downstream_frameshift", "MODIFIER"),
    }
    for julia_term, (so_term, impact) in cases.items():
        vcf = MINIMAL_HEADER + f"LmjF.01\t100\t.\tC\tT\t.\t.\tCANN=k0|.|.|{julia_term}|tx1|100|1\n"
        rows = [r for r in run_script(vcf) if r]
        assert len(rows) == 1, julia_term
        f = rows[0].split('\t')
        assert f[4] == impact, julia_term
        assert f[5] == so_term, julia_term


def test_cann_compound_effect_splits_into_one_row_per_effect():
    vcf = MINIMAL_HEADER + "LmjF.01\t100\t.\tC\tCT\t.\t.\tCANN=k0|ATG|M|missense&frameshift|tx1|100|1\n"
    rows = [r for r in run_script(vcf) if r]
    assert len(rows) == 2
    by_effect = {r.split('\t')[5]: r.split('\t')[4] for r in rows}
    assert by_effect == {"missense_variant": "MODERATE", "frameshift_variant": "HIGH"}
    assert all(r.split('\t')[6] == "product_call" for r in rows)


def test_ann_compound_effect_splits_with_replicated_impact():
    # SnpEff reports one impact for the whole entry; split rows inherit it.
    ann_compound = ("T|missense_variant&splice_region_variant|MODERATE|geneA|geneA_id|"
                    "transcript|tx1|protein_coding|1/3|c.100A>T|p.Ser34Cys|100/900|100/900|34/300|.|")
    vcf = MINIMAL_HEADER + f"LmjF.01\t3745\t.\tC\tT\t.\t.\tANN={ann_compound}\n"
    rows = [r for r in run_script(vcf) if r]
    assert len(rows) == 2
    effects = {r.split('\t')[5] for r in rows}
    assert effects == {"missense_variant", "splice_region_variant"}
    for r in rows:
        f = r.split('\t')
        assert f[4] == "MODERATE"          # entry-level impact replicated onto each row
        assert f[6] == "snpeff"


def test_cann_reference_r_entries_are_skipped():
    cann = "r0|ATG|M|reference|tx1|100|1,k0|ATG|L|missense|tx1|100|1"
    vcf  = MINIMAL_HEADER + f"LmjF.01\t100\t.\tC\tT\t.\t.\tCANN={cann}\n"
    rows = [r for r in run_script(vcf) if r]
    assert len(rows) == 1                 # only the k-entry
    assert rows[0].split('\t')[5] == "missense_variant"


def test_cann_het_indel_codon_maps_to_het_indel():
    # Codon carries an X — always a het indel (per the consensus X-masking rule).
    vcf = MINIMAL_HEADER + "LmjF.01\t100\t.\tC\tT\t.\t.\tCANN=k0|AXG|.|.|tx1|100|1\n"
    rows = [r for r in run_script(vcf) if r]
    assert len(rows) == 1
    f = rows[0].split('\t')
    assert f[5] == "het_indel"
    assert f[4] == "MODIFIER"


def test_cann_missing_sequence_codon_maps_to_coding_sequence_variant():
    # Codon is all-N (missing strain sequence, not a het indel): consequence unknown.
    vcf = MINIMAL_HEADER + "LmjF.01\t100\t.\tC\tT\t.\t.\tCANN=k0|NNN|.|.|tx1|100|1\n"
    rows = [r for r in run_script(vcf) if r]
    assert len(rows) == 1
    f = rows[0].split('\t')
    assert f[5] == "coding_sequence_variant"
    assert f[4] == "MODIFIER"


def test_header_column_names():
    vcf = MINIMAL_HEADER + f"LmjF.01\t100\t.\tC\tT\t.\t.\tANN={ANN_MODERATE}\n"
    header = run_script(vcf, include_header=True)[0]
    assert header.split('\t') == [
        "location", "seq_id", "allele", "transcript_id", "impact", "effect", "source"
    ]


def test_ann_and_cann_both_emitted_from_same_record():
    info = f"ANN={ANN_MODERATE};CANN=k0|ATG|L|missense|tx1|100|1"
    vcf  = MINIMAL_HEADER + f"LmjF.01\t3745\t.\tC\tT\t.\t.\t{info}\n"
    rows = [r for r in run_script(vcf) if r]
    assert len(rows) == 2
    sources = {r.split('\t')[6] for r in rows}
    assert sources == {"snpeff", "product_call"}


def test_cann_skips_non_coding_dot_entry():
    vcf = MINIMAL_HEADER + "LmjF.01\t100\t.\tC\tT\t.\t.\tCANN=.\n"
    rows = [r for r in run_script(vcf) if r]
    assert rows == []
