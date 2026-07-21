#!/bin/bash
#
# normalizeReadNames.sh - strip trailing mate suffixes from FASTQ read-ID lines,
# or detect whether a FASTQ carries any such suffix.
#
# bwa-mem2 pairs mate1/mate2 by requiring IDENTICAL read names; it only auto-strips
# a trailing /1 or /2. Reads carrying other mate suffixes therefore abort alignment
# with "[mem_sam_pe*] paired reads have different names". Two conventions show up in
# our source data and both trip it:
#   - SRA fastq-dump --readids  ->  .1 / .2   (e.g. SRR868698.4.1 / SRR868698.4.2)
#   - legacy Illumina/Casava    ->  :a / :b   (e.g. 90206:7:1:1:518:a / :b)
#
# Modes (first argument, default "strip"):
#   strip              read FASTQ on stdin, strip a trailing /1 /2, .1 .2, or :a :b
#                      from every read-ID line (every 4th line), write to stdout.
#                      Pure stream filter (no temp files) so it can run inside a
#                      process substitution feeding bwa-mem2 at zero extra I/O cost.
#   needs-normalizing  read FASTQ on stdin, exit 0 if any read-ID line ends in a
#                      suffix bwa CANNOT handle (.1 .2 :a :b -- NOT /1,/2, which bwa
#                      auto-strips), else exit 1. Prints nothing; it is a predicate.
#                      Full scan, no early exit: absence of a bad suffix can only be
#                      proven by reading the whole file (a single file may mix
#                      conventions), and an early awk exit would SIGPIPE the upstream
#                      zcat and, under `set -o pipefail`, invert the gate.
#
set -euo pipefail

mode="${1:-strip}"

case "$mode" in
    strip)
        exec awk 'NR % 4 == 1 { sub(/[/.:][12ab]$/, "") } 1'
        ;;
    needs-normalizing)
        exec awk 'NR % 4 == 1 && /[.:][12ab]$/ { found = 1 } END { exit !found }'
        ;;
    *)
        echo "usage: $(basename "$0") [strip|needs-normalizing]" >&2
        exit 2
        ;;
esac
