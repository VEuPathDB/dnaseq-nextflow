"""Tests for bin/parseSnpEffAnnotations.py"""
import gzip
import os
import subprocess
import tempfile

import pytest

SCRIPT = os.path.join(os.path.dirname(__file__), "../../bin/parseSnpEffAnnotations.py")


def run_script(vcf_content: str, compressed: bool = True) -> list[str]:
    """Write VCF content to a temp file, run script, return output lines (no header)."""
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
        return lines[1:]  # skip header


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
