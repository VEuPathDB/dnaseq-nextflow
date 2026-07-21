#!/usr/bin/env bash
# Unit tests for bin/normalizeReadNames.sh: the `strip` transform and the
# `needs-normalizing` predicate. No -e: `needs-normalizing` uses exit codes as
# data, and command-substitution pipelines below are inspected, not asserted.
set -uo pipefail

SCRIPT="$(cd "$(dirname "$0")/../.." && pwd)/bin/normalizeReadNames.sh"
fail() { echo "FAIL: $1"; exit 1; }

# a FASTQ record is 4 lines: @id, seq, +, qual. build one with a given id.
rec() { printf '%s\nACGT\n+\nFFFF\n' "$1"; }
# variant whose QUALITY line ends in ":1" (a suffix-looking string on a non-ID line)
recq() { printf '%s\nACGT\n+\nFF:1\n' "$1"; }

# --- strip mode --------------------------------------------------------------
out="$(rec '@read/1' | "$SCRIPT" strip)";        [ "${out%%$'\n'*}" = "@read" ]          || fail "strip /1: '${out%%$'\n'*}'"
out="$(rec '@SRR868698.4.1' | "$SCRIPT")";       [ "${out%%$'\n'*}" = "@SRR868698.4" ]   || fail "strip .1 (default mode): '${out%%$'\n'*}'"
out="$(rec '@90206:7:1:1:518:a' | "$SCRIPT" strip)"; [ "${out%%$'\n'*}" = "@90206:7:1:1:518" ] || fail "strip :a: '${out%%$'\n'*}'"
out="$(rec '@keep/1extra' | "$SCRIPT" strip)";   [ "${out%%$'\n'*}" = "@keep/1extra" ]   || fail "strip keep/1extra: '${out%%$'\n'*}'"
out="$(rec '@plain.name' | "$SCRIPT" strip)";    [ "${out%%$'\n'*}" = "@plain.name" ]    || fail "strip plain.name: '${out%%$'\n'*}'"

# only every 4th line is touched: id stripped, qual ":1" left intact
out="$(recq '@read.1' | "$SCRIPT" strip)"
[ "${out%%$'\n'*}" = "@read" ]  || fail "strip recq id: '${out%%$'\n'*}'"
[ "${out##*$'\n'}" = "FF:1" ]   || fail "strip touched a non-ID line: '${out##*$'\n'}'"

# --- needs-normalizing mode (exit 0 = needs, exit 1 = conforming) ------------
if   { rec '@read/1'; rec '@read/2'; } | "$SCRIPT" needs-normalizing; then fail "/1,/2 should be conforming (exit 1)"; fi
if   { rec '@read';   rec '@read';   } | "$SCRIPT" needs-normalizing; then fail "identical names should be conforming (exit 1)"; fi
if ! { rec '@x.1';    rec '@x.2';    } | "$SCRIPT" needs-normalizing; then fail ".1,.2 should need normalizing (exit 0)"; fi
if ! { rec '@x:a';    rec '@x:b';    } | "$SCRIPT" needs-normalizing; then fail ":a,:b should need normalizing (exit 0)"; fi

# mixed file: conforming identical-name reads (end in a coordinate, not a mate
# suffix) followed by one :a/:b pair -> whole sample needs normalizing (exit 0)
if ! { rec '@090126:3:1:1:1092'; rec '@090126:3:1:1:1092'; \
       rec '@90206:7:1:1:518:a'; rec '@90206:7:1:1:518:b'; } \
     | "$SCRIPT" needs-normalizing; then fail "mixed file should need normalizing (exit 0)"; fi

# --- usage -------------------------------------------------------------------
rc=0; "$SCRIPT" bogus </dev/null 2>/dev/null || rc=$?
[ "$rc" -eq 2 ] || fail "unknown mode should exit 2, got $rc"

echo "PASS"
