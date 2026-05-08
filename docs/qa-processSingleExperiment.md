# QA Guide: `processSingleExperiment` Output

**Purpose:** Repeatable checklist for validating the three primary outputs of a
`processSingleExperiment` run — the per-sample VCF, the consensus FASTA, and
`indels.tsv` — with focus on correctness and cross-artifact coordinate
consistency.

**Outputs under test:**

| Artifact | Path pattern |
|---|---|
| VCF | `$outputDir/<sample>.vcf.gz` |
| Consensus FASTA | `$outputDir/<sample>_consensus.fa.gz` |
| Indels table | `$outputDir/indels.tsv` |

---

## Setup

```bash
OUTPUT=/path/to/output
SAMPLE=LV39cl5           # repeat for each sample
REF=/path/to/genome.fa   # must have .fai

# Decompress and index the consensus FASTA so samtools faidx works
bgzip -d -c $OUTPUT/${SAMPLE}_consensus.fa.gz > /tmp/${SAMPLE}_consensus.fa
samtools faidx /tmp/${SAMPLE}_consensus.fa

# Helper — print a region as a single line
faidx1() { samtools faidx $1 $2 | tail -1 | tr -d '\n'; }
```

---

## Check 1: Basic Counts

Sanity-check that all three artifacts have data.

```bash
# VCF record count
zcat $OUTPUT/${SAMPLE}.vcf.gz | grep -vc "^#"

# indels.tsv lines for this sample
grep -c "^${SAMPLE}" $OUTPUT/indels.tsv

# Consensus FASTA sequence length vs reference
grep -v ">" /tmp/${SAMPLE}_consensus.fa | tr -d '\n' | wc -c
samtools faidx $REF $(grep ">" /tmp/${SAMPLE}_consensus.fa | head -1 | tr -d '>') | tail -1 | tr -d '\n' | wc -c
```

The consensus length must equal `ref_length + net_indel_offset` exactly, where
`net_indel_offset` is the sum of all `indels.tsv` entries for this sample.
No exceptions — every record in `indels.tsv` must be reflected in the consensus.

---

## Check 2: Hom SNPs → Correct ALT Base

Pick a handful of hom SNP positions and confirm the ALT base appears in the
consensus at the correct coordinate-shifted position.

```bash
# List first 10 hom SNPs
zcat $OUTPUT/${SAMPLE}.vcf.gz | grep -v "^#" \
  | awk 'substr($10,1,3)=="1/1" && $8~/TYPE=snp/' \
  | head -10 | awk '{print $1,$2,$4,$5}'

# For each position POS, compute cumulative offset (see Check 5 helper), then:
#   CONS_POS = POS + offset
# Confirm: faidx1 $REF LmjF.01:POS-POS  gives REF base
#           faidx1 /tmp/${SAMPLE}_consensus.fa LmjF.01:CONS_POS-CONS_POS  gives ALT base
```

---

## Check 3: Het SNPs → IUPAC Codes

Confirm that heterozygous SNPs are represented as IUPAC ambiguity codes, not
silently dropped or replaced with N.

```bash
# Count IUPAC codes present (excluding N)
grep -v ">" /tmp/${SAMPLE}_consensus.fa \
  | grep -io "[RYSWKMBDHV]" | sort | uniq -c | sort -rn

# Pick a 0/1 SNP and verify
zcat $OUTPUT/${SAMPLE}.vcf.gz | grep -v "^#" \
  | awk 'substr($10,1,3)=="0/1"' | grep "TYPE=snp" | head -5 \
  | awk '{print $1,$2,$4,$5}'

# For A→G het: expect R.  For C→T het: expect Y.  For G→A het: expect R.
# Reference: frozenset({A,G})=R, {C,T}=Y, {A,C}=M, {G,T}=K, {A,T}=W, {C,G}=S
```

---

## Check 4: Het Indels → REF Emitted, Not in indels.tsv

Rule: `0/1` het indels must emit REF (no coordinate shift). `1/2` het indels
must emit `X × len(REF)`.

```bash
# List het indel positions in VCF
HET_INDELS=$(zcat $OUTPUT/${SAMPLE}.vcf.gz | grep -v "^#" \
  | awk 'substr($10,1,3)=="0/1"' \
  | grep "TYPE=ins\|TYPE=del\|TYPE=complex" \
  | awk '{print $2}')

# Confirm none appear in indels.tsv
for POS in $HET_INDELS; do
    if grep -q "^${SAMPLE}	.*	${POS}	" $OUTPUT/indels.tsv; then
        echo "FAIL: het indel at $POS found in indels.tsv"
    fi
done
echo "done"

# Spot-check a few: at CONS_POS (with offset applied) expect REF bases, not X or N
# Example for het del CA→C at pos P (0/1):
#   faidx1 $REF LmjF.01:P-$((P+1))           # should show CA
#   faidx1 /tmp/${SAMPLE}_consensus.fa LmjF.01:CONS_P-$((CONS_P+1))  # should also show CA
```

---

## Check 5: Coordinate Consistency — indels.tsv vs Consensus FASTA

The cumulative offset at any reference position `X` is the sum of all
`indels.tsv` entries (for this sample) at positions `< X`. The consensus FASTA
coordinate for reference position `X` is `X + offset`.

```bash
# Compute offset at a given reference position REF_POS
cumulative_offset() {
    local ref_pos=$1
    grep "^${SAMPLE}	LmjF.01" $OUTPUT/indels.tsv \
      | awk -v p=$ref_pos '$3 < p {sum += $4} END {print sum+0}'
}

# Verify a hom indel: pick a deletion from indels.tsv, confirm the bases it
# deletes are absent from the consensus and the next base aligns correctly.
# Example: hom del CGA→C at pos P.  offset before P = O.  cons pos for P = P+O.
#   REF  at P:(P+2)  → should show CGA
#   CONS at (P+O):(P+O)  → should show C (anchor only, GA gone)
#   CONS at (P+O+1)  → should be ref base at P+3

# Accumulate and print offset at each indel position to spot jumps
grep "^${SAMPLE}	LmjF.01" $OUTPUT/indels.tsv \
  | awk '{sum += $4; print $3, $4, "cumsum:", sum}'
```

---

## Check 6: Complex Variants — Atomized Published VCF vs Single Consensus Record

FreeBayes emits a single TYPE=complex record for sites with both a SNP and an
indel component. The pipeline branches these into two paths:

- **`consensus_input.vcf.gz`** — retains the raw complex record as-is; the
  consensus builder applies the full allele in one step.
- **`<sample>.vcf.gz`** (published) — atomized: the complex record is split
  into separate SNP and indel rows at the same position.

This means the published VCF will have **2 rows** at a complex site while the
consensus was built from **1 row**. This is expected and correct — the indel
offset applied to the consensus and recorded in `indels.tsv` are consistent.

**How to verify a complex site:**

```bash
WD=<filterAndSplitVcf work dir>

# consensus_input: should show 1 complex record
bcftools view $WD/consensus_input.vcf.gz ${CHROM}:${POS}-${POS} | grep -v '^#'

# published VCF: should show 2 atomized records at the same position
bcftools view $WD/${SAMPLE}.vcf.gz ${CHROM}:${POS}-${POS} | grep -v '^#'
```

**Co-located hom records in the published VCF are not a bug** — they are the
expected result of atomization. The Check 1 length invariant
(`consensus = ref + net_indel_offset`) must still hold exactly.

---

## Check 7: Hom Large Deletions — Sequence Gap in Consensus

For large deletions confirm the deleted bases are absent and the anchor plus
downstream sequence align to the reference.

```bash
# Find large hom deletions (> 5 bp)
zcat $OUTPUT/${SAMPLE}.vcf.gz | grep -v "^#" \
  | awk 'substr($10,1,3)=="1/1" && length($4) > 5 && length($5)==1' \
  | awk '{print $1,$2,$4,$5,length($4)-1,"bp deleted"}'

# Spot-check: for del REF=NNNNN ALT=N at pos P, offset O = cumulative before P
# faidx1 $REF LmjF.01:(P-2)-(P+length(REF)+2)   # shows full context in ref
# faidx1 /tmp/${SAMPLE}_consensus.fa LmjF.01:(P+O-2)-(P+O+2)  # deleted bases gone, anchor present
```

---

## Check 8: Coverage Masking — Uncovered Regions → N

Non-variant positions with zero coverage (gaps between BED intervals) should
be `N` in the consensus. The BED file only lists covered windows; positions not
covered by any interval are gaps.

```bash
# Find gaps between BED intervals (positions not covered by any row)
zcat $OUTPUT/${SAMPLE}.coverage.bed.gz \
  | awk '{print $2,$3}' | sort -n \
  | awk 'NR==1{prev=$2;next} $1>prev{print "gap:", prev, $1, "len:", $1-prev; prev=$2} {prev=$2}' \
  | head -5

# For a gap starting at GAP_START (0-based), compute the consensus position
# accounting for cumulative indel offset, then confirm it is N:
OFFSET=$(grep "^${SAMPLE}	LmjF.01" $OUTPUT/indels.tsv \
  | awk -v p=$((GAP_START+1)) '$3 < p {sum+=$4} END {print sum+0}')
CONS_POS=$(( GAP_START + 1 + OFFSET ))
faidx1 /tmp/${SAMPLE}_consensus.fa LmjF.01:${CONS_POS}-${CONS_POS}
# expect: N
```

**Important:** the BED coverage check is only applied to non-variant gaps.
FreeBayes enforces `--min-coverage` before emitting any record, so VCF
positions are never re-checked against the BED. A hom indel at a position that
bedtools reports as zero-coverage (because all reads carry the deletion) is
still applied correctly — the called allele appears in the consensus and the
offset is recorded in `indels.tsv`.

---

## Check 9: Seidman751 / Multi-sample Consistency

When multiple samples are present, confirm `indels.tsv` has the correct sample
column and that coordinate checks for each sample use that sample's own offsets.

```bash
# Confirm sample names present
cut -f1 $OUTPUT/indels.tsv | sort -u

# Re-run Checks 1–8 substituting SAMPLE=Seidman751 (or each additional sample)
```

---

## Check 10: Het SNPs — Exhaustive IUPAC Verification

Every `0/1` SNP in the published VCF must produce the correct IUPAC ambiguity code
at the coordinate-shifted consensus position. This is an exhaustive check, not a
spot-check.

IUPAC mapping used: `{A,G}=R  {C,T}=Y  {G,C}=S  {A,T}=W  {G,T}=K  {A,C}=M`

```bash
OUTPUT=/path/to/output
REF=/path/to/genome.fa
SAMPLE=LV39cl5

declare -A IUPAC
IUPAC[AG]=R; IUPAC[GA]=R
IUPAC[CT]=Y; IUPAC[TC]=Y
IUPAC[GC]=S; IUPAC[CG]=S
IUPAC[AT]=W; IUPAC[TA]=W
IUPAC[GT]=K; IUPAC[TG]=K
IUPAC[AC]=M; IUPAC[CA]=M

PASS=0; FAIL=0
while IFS=$'\t' read -r CHROM POS ID REF_BASE ALT_BASE REST; do
    OFFSET=$(grep "^${SAMPLE}	${CHROM}" $OUTPUT/indels.tsv \
        | awk -v p="$POS" '$3 < p {sum += $4} END {print sum+0}')
    CONS_POS=$((POS + OFFSET))
    CONS_BASE=$(samtools faidx /tmp/${SAMPLE}_consensus.fa ${CHROM}:${CONS_POS}-${CONS_POS} | tail -1 | tr 'a-z' 'A-Z')
    REF_U=$(echo "$REF_BASE" | tr 'a-z' 'A-Z')
    ALT_U=$(echo "$ALT_BASE" | tr 'a-z' 'A-Z')
    EXPECT="${IUPAC[${REF_U}${ALT_U}]}"
    if [ "$CONS_BASE" = "$EXPECT" ]; then PASS=$((PASS+1))
    else echo "FAIL: POS=$POS REF=$REF_U ALT=$ALT_U expect=$EXPECT cons=$CONS_BASE"; FAIL=$((FAIL+1))
    fi
done < <(zcat $OUTPUT/${SAMPLE}.vcf.gz | grep -v "^#" \
    | awk 'substr($10,1,3)=="0/1" && $8~/TYPE=snp/' | cut -f1-5,7-)

echo "PASS: $PASS  FAIL: $FAIL"
# Repeat for each SAMPLE
```

All `FAIL` lines indicate a bug. `PASS` count should equal the total `0/1 TYPE=snp`
record count in the VCF.

**Note:** Atomized complex variants carry `TYPE=complex` in INFO even for their
single-base SNP row — those are caught by Check 11, not here.

---

## Check 11: No Duplicate SNP Rows — Multi-row Sites Are SNP + Indel Only

A SNP call (single-base REF → single-base ALT) must never appear on more than one
row at the same position. The only valid reason for two rows at the same position
is an atomized complex variant: one SNP-like row (`len(REF)==1 && len(ALT)==1`) and
one indel-like row (`len(REF) != len(ALT)`).

**Why `TYPE=` is not the right filter:** atomized complex variants retain
`TYPE=complex` in INFO for both the SNP and indel rows. Length-based discrimination
is the correct approach.

```bash
OUTPUT=/path/to/output
SAMPLE=LV39cl5

# Find positions with more than one row
MULTI=$(zcat $OUTPUT/${SAMPLE}.vcf.gz | grep -v "^#" \
    | awk '{print $1"_"$2}' | sort | uniq -d)
echo "Positions with >1 row: $(echo "$MULTI" | grep -c .)"

COUNT_GOOD=0; COUNT_BAD=0
for SITE in $MULTI; do
    CHROM=$(echo $SITE | cut -d_ -f1)
    POS=$(echo $SITE | cut -d_ -f2)
    SNP_LIKE=$(zcat $OUTPUT/${SAMPLE}.vcf.gz | grep -v "^#" \
        | awk -v c=$CHROM -v p=$POS '$1==c && $2==p && length($4)==1 && length($5)==1' | wc -l)
    TOTAL=$(zcat $OUTPUT/${SAMPLE}.vcf.gz | grep -v "^#" \
        | awk -v c=$CHROM -v p=$POS '$1==c && $2==p' | wc -l)
    if [ "$SNP_LIKE" -le 1 ] && [ "$TOTAL" -eq 2 ]; then
        COUNT_GOOD=$((COUNT_GOOD+1))
    else
        COUNT_BAD=$((COUNT_BAD+1))
        echo "FAIL $SITE: total=$TOTAL snp-like=$SNP_LIKE"
        zcat $OUTPUT/${SAMPLE}.vcf.gz | grep -v "^#" \
            | awk -v c=$CHROM -v p=$POS '$1==c && $2==p {print "  ref="$4,"alt="$5}'
    fi
done
echo "Good (1 SNP-like + 1 indel-like): $COUNT_GOOD  Bad: $COUNT_BAD"
# Repeat for each SAMPLE
```

Any `FAIL` line means a SNP was emitted on two separate rows at the same position,
which would cause double-counting in downstream merge.

---

## Quick-pass Summary Template

```
Sample: _______________   Run date: _______________

[ ] Check 1: VCF rows > 0, indels.tsv rows > 0, consensus length reasonable
[ ] Check 2: Hom SNP base correct at shifted consensus position
[ ] Check 3: IUPAC codes present (R/Y/S/K/W/M counts > 0)
[ ] Check 4: Het indels absent from indels.tsv; REF sequence at het indel sites
[ ] Check 5: Cumulative offset consistent between indels.tsv and consensus coordinates
[ ] Check 6: Co-located complex SNP+indel positions identified and phantom offset quantified
[ ] Check 7: Large deletions absent from consensus (spot-checked ≥ 1)
[ ] Check 8: Zero-coverage regions are N in consensus
[ ] Check 9: Multi-sample: all samples present in indels.tsv
[ ] Check 10: All het SNPs (0/1 TYPE=snp) map to correct IUPAC code — PASS count = total het SNP count
[ ] Check 11: No duplicate SNP rows; all multi-row positions are exactly 1 SNP-like + 1 indel-like

Notes / failures:
```
