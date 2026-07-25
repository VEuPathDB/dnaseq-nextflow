#!/usr/bin/env bash
# Unit tests for bin/normalizePairedReadNames.sh.
#
# The contract under test is the invariant bwa-mem2 (and samtools fixmate) needs:
# mate1[i] and mate2[i] must carry IDENTICAL read names. The normalizer asserts
# that per pair, repairing it only by removing a trailing mate suffix, and aborts
# rather than emit silently mispaired output.
#
# No -e: several cases assert on non-zero exit codes.
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SCRIPT="$ROOT/bin/normalizePairedReadNames.sh"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

fail() { echo "FAIL: $1"; exit 1; }

# a FASTQ record is 4 lines: @id, seq, +, qual
rec()  { printf '%s\nACGT\n+\nFFFF\n' "$1"; }
# variant whose QUALITY line ends in ":1" -- a suffix-looking string off the ID line
recq() { printf '%s\nACGT\n+\nFF:1\n' "$1"; }

# run the normalizer on two stdin-built mates; sets rc, OUT1, OUT2, ERR
run() {
  printf '%s' "$1" > "$TMP/in1"
  printf '%s' "$2" > "$TMP/in2"
  rc=0
  "$SCRIPT" "$TMP/in1" "$TMP/in2" "$TMP/out1" "$TMP/out2" 2>"$TMP/err" || rc=$?
  OUT1="$(cat "$TMP/out1" 2>/dev/null)"
  OUT2="$(cat "$TMP/out2" 2>/dev/null)"
  ERR="$(cat "$TMP/err")"
}

# --- SRA fastq-dump --readids: mate suffix followed by a comment ---------------
# This is the case the previous end-of-line-anchored implementation missed: the
# line ends in "length=40", so the suffix is not at end of line.
run "$(rec '@SRR868698.1.1 WIGTC-HISEQ:2:1101:1126:2062 length=40')" \
    "$(rec '@SRR868698.1.2 WIGTC-HISEQ:2:1101:1126:2062 length=53')"
[ "$rc" -eq 0 ] || fail "SRA readids+comment: exit $rc ($ERR)"
[ "${OUT1%%$'\n'*}" = '@SRR868698.1 WIGTC-HISEQ:2:1101:1126:2062 length=40' ] \
  || fail "SRA readids+comment mate1: '${OUT1%%$'\n'*}'"
[ "${OUT2%%$'\n'*}" = '@SRR868698.1 WIGTC-HISEQ:2:1101:1126:2062 length=53' ] \
  || fail "SRA readids+comment mate2: '${OUT2%%$'\n'*}'"

# --- legacy Illumina :a/:b, with and without a comment ------------------------
run "$(rec '@90206:7:1:1:518:a')" "$(rec '@90206:7:1:1:518:b')"
[ "$rc" -eq 0 ] && [ "${OUT1%%$'\n'*}" = '@90206:7:1:1:518' ] \
  || fail ":a/:b bare: rc=$rc '${OUT1%%$'\n'*}' ($ERR)"

run "$(rec '@90206:7:1:1:518:a some comment')" "$(rec '@90206:7:1:1:518:b some comment')"
[ "$rc" -eq 0 ] && [ "${OUT1%%$'\n'*}" = '@90206:7:1:1:518 some comment' ] \
  || fail ":a/:b with comment: rc=$rc '${OUT1%%$'\n'*}' ($ERR)"

# --- /1,/2: bwa strips these itself, but fixmate still needs equal names ------
run "$(rec '@read/1 c')" "$(rec '@read/2 c')"
[ "$rc" -eq 0 ] && [ "${OUT1%%$'\n'*}" = '@read c' ] \
  || fail "/1,/2: rc=$rc '${OUT1%%$'\n'*}' ($ERR)"

# --- names ALREADY identical are passed through byte-for-byte ------------------
# Critical no-false-positive case: plain fastq-dump (no --readids) numbers spots,
# so mate1 and mate2 read 1 are both "@SRR868698.1". That trailing ".1" is a spot
# index, NOT a mate suffix, and must survive.
same="$(rec '@SRR868698.1 WIGTC-HISEQ:2:1101:20786:2000 length=40')"
run "$same" "$same"
[ "$rc" -eq 0 ] || fail "identical spot-id names: exit $rc ($ERR)"
[ "$OUT1" = "${same%$'\n'}" ] || fail "identical spot-id names were altered: '$OUT1'"
[ "$OUT2" = "${same%$'\n'}" ] || fail "identical spot-id names were altered (mate2): '$OUT2'"

run "$(rec '@read')" "$(rec '@read')"
[ "$rc" -eq 0 ] && [ "${OUT1%%$'\n'*}" = '@read' ] \
  || fail "identical bare names: rc=$rc '${OUT1%%$'\n'*}' ($ERR)"

# --- only the ID line is rewritten -------------------------------------------
run "$(recq '@x.1')" "$(recq '@x.2')"
[ "$rc" -eq 0 ] || fail "recq: exit $rc ($ERR)"
[ "${OUT1%%$'\n'*}" = '@x' ]  || fail "recq id: '${OUT1%%$'\n'*}'"
[ "${OUT1##*$'\n'}" = 'FF:1' ] || fail "normalizer touched a non-ID line: '${OUT1##*$'\n'}'"

# --- mixed conventions inside one file are handled per pair -------------------
run "$(rec '@090126:3:1:1:1092'; rec '@x.1'; rec '@y/1')" \
    "$(rec '@090126:3:1:1:1092'; rec '@x.2'; rec '@y/2')"
[ "$rc" -eq 0 ] || fail "mixed conventions: exit $rc ($ERR)"
ids="$(awk 'NR%4==1{print $1}' "$TMP/out1" | tr '\n' ',')"
[ "$ids" = '@090126:3:1:1:1092,@x,@y,' ] || fail "mixed conventions ids: '$ids'"

# --- irreconcilable mismatch aborts loud, naming both reads ------------------
run "$(rec '@readA.1')" "$(rec '@readB.2')"
[ "$rc" -ne 0 ] || fail "desynced pair should abort, got exit 0"
case "$ERR" in *readA.1*readB.2*) ;; *) fail "abort message lacks both names: '$ERR'" ;; esac

# a suffix mismatch that is not a mate suffix is a mismatch, not a repair
run "$(rec '@read.3')" "$(rec '@read.4')"
[ "$rc" -ne 0 ] || fail ".3/.4 is not a mate suffix; should abort, got exit 0"

# stripping must not be allowed to produce an empty read name
run "$(rec '@.1')" "$(rec '@.2')"
[ "$rc" -ne 0 ] || fail "empty-after-strip name should abort, got exit 0"

# --- truncation / desync in record count aborts ------------------------------
run "$(rec '@x.1'; rec '@z.1')" "$(rec '@x.2')"
[ "$rc" -ne 0 ] || fail "mate2 ending early should abort, got exit 0"

run "$(rec '@x.1')" "$(rec '@x.2'; rec '@z.2')"
[ "$rc" -ne 0 ] || fail "mate1 ending early should abort, got exit 0"

# a partial 4-line record is corruption, not a valid stream
run "$(printf '@x.1\nACGT\n')" "$(printf '@x.2\nACGT\n')"
[ "$rc" -ne 0 ] || fail "truncated record should abort, got exit 0"

# --- record-frame corruption is caught --------------------------------------
run "$(printf '@x.1\nACGT\n+\nFFFF\nACGT\n+\nFFFF\n@y.1\n')" \
    "$(printf '@x.2\nACGT\n+\nFFFF\nACGT\n+\nFFFF\n@y.2\n')"
[ "$rc" -ne 0 ] || fail "ID line not starting with @ should abort, got exit 0"

# --- usage ------------------------------------------------------------------
rc=0; "$SCRIPT" only-one-arg >/dev/null 2>&1 || rc=$?
[ "$rc" -eq 2 ] || fail "wrong argument count should exit 2, got $rc"

echo "PASS"
