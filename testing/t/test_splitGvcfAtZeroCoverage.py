import os
import sys
import tempfile
import textwrap

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../bin'))

from splitGvcfAtZeroCoverage import load_zero_cov_bed, covered_intervals, process_gvcf


def test_load_zero_cov_bed_parses_intervals():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write('chr1\t100\t200\n')
        f.write('chr1\t500\t600\n')
        f.write('chr2\t0\t50\n')
        fname = f.name
    try:
        result = load_zero_cov_bed(fname)
        assert result['chr1'] == [(100, 200), (500, 600)]
        assert result['chr2'] == [(0, 50)]
    finally:
        os.unlink(fname)


def test_load_zero_cov_bed_empty_file():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        fname = f.name
    try:
        result = load_zero_cov_bed(fname)
        assert result == {}
    finally:
        os.unlink(fname)


def test_process_gvcf_splits_ref_block_at_zero_coverage():
    """A reference block spanning a zero-coverage gap should be split into two sub-blocks."""
    gvcf_content = textwrap.dedent("""\
        ##fileformat=VCFv4.2
        ##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
        chr1\t1\t.\tA\t<*>\t.\t.\tEND=1000\tGT:DP\t0/0:20\t0/0:15
    """)
    bed_content = 'chr1\t499\t501\n'  # zero-coverage gap at positions 500-501

    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as gvcf_f, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as bed_f, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as out_f:
        gvcf_f.write(gvcf_content)
        bed_f.write(bed_content)
        gvcf_path = gvcf_f.name
        bed_path = bed_f.name
        out_path = out_f.name

    try:
        process_gvcf(gvcf_path, bed_path, out_path)
        with open(out_path) as f:
            lines = [l for l in f if not l.startswith('#')]
        assert len(lines) == 2, f"Expected 2 sub-blocks, got {len(lines)}: {lines}"
        pos1 = int(lines[0].split('\t')[1])
        pos2 = int(lines[1].split('\t')[1])
        assert pos1 == 1
        assert pos2 == 502  # first base after zero-cov gap (501 0-based = 502 1-based)
    finally:
        for p in [gvcf_path, bed_path, out_path]:
            os.unlink(p)


def test_process_gvcf_drops_fully_zero_block():
    """A reference block entirely within a zero-coverage region should be dropped."""
    gvcf_content = textwrap.dedent("""\
        ##fileformat=VCFv4.2
        ##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
        chr1\t100\t.\tA\t<*>\t.\t.\tEND=200\tGT:DP\t0/0:0
    """)
    bed_content = 'chr1\t99\t201\n'  # covers entire block (0-based 99-201 = 1-based 100-201)

    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as gvcf_f, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as bed_f, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as out_f:
        gvcf_f.write(gvcf_content)
        bed_f.write(bed_content)
        gvcf_path = gvcf_f.name
        bed_path = bed_f.name
        out_path = out_f.name

    try:
        process_gvcf(gvcf_path, bed_path, out_path)
        with open(out_path) as f:
            lines = [l for l in f if not l.startswith('#')]
        assert len(lines) == 0, f"Expected block to be dropped, got: {lines}"
    finally:
        for p in [gvcf_path, bed_path, out_path]:
            os.unlink(p)


def test_process_gvcf_passes_variant_records_unchanged():
    """Non-reference-block records should pass through untouched."""
    gvcf_content = textwrap.dedent("""\
        ##fileformat=VCFv4.2
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
        chr1\t100\t.\tA\tT\t50\t.\t.\tGT:DP\t1/1:30\t0/1:25
    """)
    bed_content = 'chr1\t99\t101\n'

    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as gvcf_f, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as bed_f, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as out_f:
        gvcf_f.write(gvcf_content)
        bed_f.write(bed_content)
        gvcf_path = gvcf_f.name
        bed_path = bed_f.name
        out_path = out_f.name

    try:
        process_gvcf(gvcf_path, bed_path, out_path)
        with open(out_path) as f:
            lines = [l for l in f if not l.startswith('#')]
        assert len(lines) == 1
        assert lines[0].split('\t')[1] == '100'
    finally:
        for p in [gvcf_path, bed_path, out_path]:
            os.unlink(p)
