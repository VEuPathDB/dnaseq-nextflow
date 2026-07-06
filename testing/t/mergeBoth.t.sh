#!/usr/bin/env bash
# Characterization: bcftools merge -m both keeps SNP and indel alleles on
# SEPARATE records at a shared start position (never a mixed multiallelic row).
set -euo pipefail
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT; cd "$tmp"

cat > hdr.txt <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
EOF
mk() { { cat hdr.txt; printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$1";
         printf 'chr1\t100\t.\t%s\t%s\t.\t.\t.\tGT\t1/1\n' "$2" "$3"; } > "$1.vcf"
       bgzip -f "$1.vcf"; tabix -f -p vcf "$1.vcf.gz"; }
mk SA A G       # SNP
mk SB A AT      # insertion

out="$(bcftools merge -m both SA.vcf.gz SB.vcf.gz 2>/dev/null | grep -vc '^#')"
[ "$out" = "2" ] || { echo "FAIL: expected 2 records, got $out"; exit 1; }

allrows="$(bcftools merge -m all SA.vcf.gz SB.vcf.gz 2>/dev/null | grep -vc '^#')"
[ "$allrows" = "1" ] || { echo "FAIL: expected -m all to fuse into 1 row, got $allrows"; exit 1; }

echo "PASS"
