#!/usr/bin/env python3
"""Compare two mergeExperiments output directories for equivalence modulo
strain-id permutation.

Strain ids are assigned in whatever order the collected per-strain files were
staged, so two runs over identical inputs can produce identical data under
different numbering. This canonicalizes every id-bearing field to sample NAMES
before comparing, and sorts rows (row order within a file is not guaranteed
either).

Usage: compareMergeOutputs.py DIR_A DIR_B
Exit 0 if equivalent, 1 if not, 2 on a usage/IO problem.
"""

import subprocess
import sys
from pathlib import Path

# .dat files with no id-bearing column: compared as-is (sorted).
PLAIN_DATS = ["variationFeature.dat", "snpeff.dat"]

# .dat files with a column holding a "{1,3}"-style strain-id set.
ID_SET_DATS = {
    "allele.dat": "strain_ids",
    "transcript_product.dat": "downstream_of_frameshift_strain_ids",
}

problems = []


def problem(msg):
    problems.append(msg)


def read_id_map(path):
    """Read a two-column id<TAB>name file into {id: name}, skipping a header
    row if present. Used for both sample.dat and strainIdToName.dat."""
    mapping = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            problem(f"{path}: cannot parse line: {line!r}")
            continue
        sid, name = parts[0], parts[1]
        if sid == "strain_id":  # sample.dat header
            continue
        mapping[sid] = name
    return mapping


def canon_id_set(field, id_map, where):
    """'{1,3}' -> '{Friedlin,LV39}' (name-sorted). Empty stays empty."""
    raw = field.strip()
    if not raw or raw == ".":
        return raw
    if not (raw.startswith("{") and raw.endswith("}")):
        problem(f"{where}: expected a {{...}} id set, got {field!r}")
        return raw
    ids = [i for i in raw[1:-1].split(",") if i]
    names = []
    for i in ids:
        if i not in id_map:
            problem(f"{where}: strain id {i} not in sample.dat")
            names.append(f"<unknown:{i}>")
        else:
            names.append(id_map[i])
    return "{" + ",".join(sorted(names)) + "}"


def load_dat(path, id_col=None, id_map=None):
    """Return (header, sorted canonicalized data rows)."""
    lines = path.read_text().splitlines()
    if not lines:
        return "", []
    header, rows = lines[0], lines[1:]
    if id_col is None:
        return header, sorted(rows)
    cols = header.lstrip("#").split("\t")
    if id_col not in cols:
        problem(f"{path}: no {id_col} column in header")
        return header, sorted(rows)
    idx = cols.index(id_col)
    out = []
    for row in rows:
        fields = row.split("\t")
        if idx < len(fields):
            fields[idx] = canon_id_set(fields[idx], id_map, f"{path.name}")
        out.append("\t".join(fields))
    return header, sorted(out)


def compare_rows(name, a, b):
    ah, arows = a
    bh, brows = b
    if ah != bh:
        problem(f"{name}: header differs\n  A: {ah}\n  B: {bh}")
    if arows == brows:
        return
    only_a = [r for r in arows if r not in set(brows)]
    only_b = [r for r in brows if r not in set(arows)]
    problem(
        f"{name}: {len(only_a)} row(s) only in A, {len(only_b)} row(s) only in B\n"
        + "".join(f"  A only: {r}\n" for r in only_a[:5])
        + "".join(f"  B only: {r}\n" for r in only_b[:5])
    )


def compare_sample_dat(a_dir, b_dir):
    """Sample identity is the set of names; the id assignment is what we are
    deliberately ignoring."""
    a = set(read_id_map(a_dir / "sample.dat").values())
    b = set(read_id_map(b_dir / "sample.dat").values())
    if a != b:
        problem(f"sample.dat: sample sets differ\n  only A: {sorted(a - b)}\n  only B: {sorted(b - a)}")


def compare_hsss(a_dir, b_dir):
    """hsss dirs hold per-strain files NAMED by a strain index that uses its own
    numbering (reference first). Map filenames to names via each dir's own
    strainIdToName.dat, then compare contents pairwise by name."""
    for freq in ("20", "40", "60", "80"):
        d = f"hsss_readFreq{freq}"
        a, b = a_dir / d, b_dir / d
        if not a.is_dir() or not b.is_dir():
            problem(f"{d}: missing in {'A' if not a.is_dir() else 'B'}")
            continue
        for meta in ("contigIdToSourceId.dat", "referenceGenome.dat"):
            if (a / meta).read_bytes() != (b / meta).read_bytes():
                problem(f"{d}/{meta}: differs")
        a_map = read_id_map(a / "strainIdToName.dat")
        b_map = read_id_map(b / "strainIdToName.dat")
        if set(a_map.values()) != set(b_map.values()):
            problem(f"{d}/strainIdToName.dat: strain sets differ")
            continue
        b_by_name = {name: sid for sid, name in b_map.items()}
        for sid, name in a_map.items():
            af, bf = a / sid, b / b_by_name[name]
            if not af.exists() or not bf.exists():
                problem(f"{d}: no per-strain file for {name} (A:{af.name} B:{bf.name})")
                continue
            if af.read_bytes() != bf.read_bytes():
                problem(f"{d}: per-strain data for {name} differs ({af.name} vs {bf.name})")


def read_vcf(path):
    """Return (sample_names, record_lines), reading the text directly rather
    than through bcftools.

    Two reasons not to use bcftools here. First, the `##` header carries
    bcftools and snpEff version and command-line stamps, so it must be skipped
    entirely -- comparing raw bytes of the compressed file is meaningless.
    Second, this pipeline's published VCF uses FORMAT tags it never declares
    (rows carry GT:DP:AD:RO:QR:AO:QA:CA:DFS while only CA and DFS are declared),
    which makes any bcftools mode that validates the header -- including the
    `view -s` sample subsetting this function used to rely on -- refuse to run.
    Reordering the sample columns ourselves sidesteps that entirely.
    """
    out = subprocess.run(["bgzip", "-d", "-c", str(path)],
                         capture_output=True, text=True, check=True).stdout
    samples, records = [], []
    for line in out.splitlines():
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            samples = line.split("\t")[9:]
            continue
        records.append(line)
    return samples, records


def canon_vcf_records(samples, records):
    """Records with sample columns permuted into name-sorted order, then sorted.
    This is what makes the comparison blind to sample column order."""
    order = sorted(range(len(samples)), key=lambda i: samples[i])
    out = []
    for line in records:
        fields = line.split("\t")
        fixed, cols = fields[:9], fields[9:]
        out.append("\t".join(fixed + [cols[i] for i in order]))
    return sorted(out)


def compare_ann_vcf(a_dir, b_dir):
    name = "merged.ann.vcf.gz"
    a, b = a_dir / name, b_dir / name
    if not a.exists() or not b.exists():
        problem(f"{name}: missing in {'A' if not a.exists() else 'B'}")
        return

    try:
        a_samples, a_records = read_vcf(a)
        b_samples, b_records = read_vcf(b)
    except subprocess.CalledProcessError as e:
        problem(f"{name}: could not read ({e.cmd[0]} exited {e.returncode}): {e.stderr.strip()}")
        return

    if sorted(a_samples) != sorted(b_samples):
        problem(f"{name}: sample sets differ\n  A: {sorted(a_samples)}\n  B: {sorted(b_samples)}")
        return
    if canon_vcf_records(a_samples, a_records) != canon_vcf_records(b_samples, b_records):
        problem(f"{name}: record body differs "
                f"({len(a_records)} records in A, {len(b_records)} in B)")


def main():
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    a_dir, b_dir = Path(sys.argv[1]), Path(sys.argv[2])
    for d in (a_dir, b_dir):
        if not d.is_dir():
            print(f"not a directory: {d}", file=sys.stderr)
            return 2

    a_map = read_id_map(a_dir / "sample.dat")
    b_map = read_id_map(b_dir / "sample.dat")

    compare_sample_dat(a_dir, b_dir)

    for name in PLAIN_DATS:
        compare_rows(name, load_dat(a_dir / name), load_dat(b_dir / name))

    for name, col in ID_SET_DATS.items():
        compare_rows(name,
                     load_dat(a_dir / name, col, a_map),
                     load_dat(b_dir / name, col, b_map))

    compare_hsss(a_dir, b_dir)
    compare_ann_vcf(a_dir, b_dir)

    if problems:
        print(f"DIFFER: {len(problems)} problem(s)\n")
        for p in problems:
            print(f"- {p}")
        return 1
    print("EQUIVALENT (modulo strain-id permutation)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
