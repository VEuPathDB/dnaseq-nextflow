#!/bin/bash
#
# normalizeReadNames.sh - strip trailing mate suffixes from FASTQ read-ID lines.
#
# bwa-mem2 pairs mate1/mate2 by requiring IDENTICAL read names; it only auto-strips
# a trailing /1 or /2. Reads carrying other mate suffixes therefore abort alignment
# with "[mem_sam_pe*] paired reads have different names". Two conventions show up in
# our source data and both trip it:
#   - SRA fastq-dump --readids  ->  .1 / .2   (e.g. SRR868698.4.1 / SRR868698.4.2)
#   - legacy Illumina/Casava    ->  :a / :b   (e.g. 90206:7:1:1:518:a / :b)
#
# This filter reads FASTQ on stdin and writes it to stdout, stripping a trailing
# /1 /2, .1 .2, or :a :b from every read-ID line (every 4th line) so mate names match.
# It is a pure stream filter (no temp files) so it can run inside a process
# substitution feeding bwa-mem2 at zero extra I/O cost.
#
set -euo pipefail
exec awk 'NR % 4 == 1 { sub(/[/.:][12ab]$/, "") } 1'
