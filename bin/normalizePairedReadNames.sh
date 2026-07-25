#!/bin/bash
#
# normalizePairedReadNames.sh IN1 IN2 OUT1 OUT2
#
# Enforce the invariant paired alignment actually needs: for every read pair,
# mate1 and mate2 must carry the SAME read name. bwa-mem2 pairs mates by name
# equality (auto-stripping only a trailing /1 or /2) and aborts with
# "[mem_sam_pe*] paired reads have different names"; samtools fixmate likewise
# groups mates by QNAME. Source data breaks the invariant in two ways:
#   - SRA fastq-dump --readids ->  .1 / .2   (@SRR868698.4.1 / @SRR868698.4.2)
#   - legacy Illumina/Casava   ->  :a / :b   (@90206:7:1:1:518:a / :b)
#
# This is a PAIRWISE property, so it is checked pairwise rather than guessed at
# from the lexical shape of one file's read names. Both mates are streamed in
# lockstep, one pass, and for each pair:
#   names equal                                  -> emit both records verbatim
#   equal after removing a trailing /1 /2 .1 .2
#     :a :b from each name                       -> emit with the suffix removed
#   anything else (or a record-count / framing
#     mismatch between the mates)                -> abort, non-zero exit
#
# Two consequences worth knowing:
#   - Names that are already identical are passed through untouched. This
#     matters: plain fastq-dump (no --readids) numbers SPOTS, so read 1 of both
#     mates is "@SRR868698.1". That trailing ".1" is not a mate suffix and a
#     purely lexical stripper would mangle it.
#   - Only the read-name field (up to the first whitespace) is considered. SRA
#     and Casava headers carry a trailing comment ("... length=40"), so the mate
#     suffix is mid-line, not at end of line.
#
# IN1/IN2 and OUT1/OUT2 may be FIFOs or /dev/fd paths, so this composes with
# process substitution and needs no temporary files. Aborting mid-stream leaves
# partially written output: callers MUST treat a non-zero exit as fatal rather
# than trusting whatever the downstream aligner made of a truncated stream.
#
set -uo pipefail

if [ "$#" -ne 4 ]; then
    echo "usage: $(basename "$0") IN1 IN2 OUT1 OUT2" >&2
    exit 2
fi

exec awk -v f1="$1" -v f2="$2" -v o1="$3" -v o2="$4" '
function die(msg) {
    printf("ERROR: normalizePairedReadNames: %s\n", msg) | "cat 1>&2"
    aborted = 1
    exit 1
}
# read the 3 non-ID lines of a record from each mate, keeping them verbatim
function readBody(pair,   i, l1, l2) {
    for (i = 1; i <= 3; i++) {
        if ((getline l1 < f1) <= 0 || (getline l2 < f2) <= 0)
            die("record " pair " is truncated (a FASTQ record is 4 lines)")
        body1[i] = l1
        body2[i] = l2
    }
}
BEGIN {
    while (1) {
        r1 = (getline h1 < f1)
        r2 = (getline h2 < f2)
        if (r1 < 0) die("cannot read " f1)
        if (r2 < 0) die("cannot read " f2)
        if (r1 == 0 && r2 == 0) break
        if (r1 == 0) die("mate 1 ended after " pairs " reads while mate 2 continues; mates are not in lockstep")
        if (r2 == 0) die("mate 2 ended after " pairs " reads while mate 1 continues; mates are not in lockstep")

        pairs++
        if (substr(h1, 1, 1) != "@" || substr(h2, 1, 1) != "@")
            die("record " pairs " does not start with @; FASTQ framing is corrupt")
        readBody(pairs)

        # split each header into the read-name field and its trailing comment
        match(h1, /^[^ \t]+/); n1 = substr(h1, 1, RLENGTH); rest1 = substr(h1, RLENGTH + 1)
        match(h2, /^[^ \t]+/); n2 = substr(h2, 1, RLENGTH); rest2 = substr(h2, RLENGTH + 1)

        if (n1 != n2) {
            s1 = n1; sub(/[/.:][12ab]$/, "", s1)
            s2 = n2; sub(/[/.:][12ab]$/, "", s2)
            if (s1 != s2)
                die("read " pairs " names differ irreconcilably: \"" substr(n1, 2) "\" vs \"" \
                    substr(n2, 2) "\"; mates are mispaired, not merely mate-suffixed")
            if (length(s1) < 2)
                die("read " pairs " names \"" substr(n1, 2) "\" / \"" substr(n2, 2) \
                    "\" are nothing but a mate suffix; stripping would leave an empty read name")
            n1 = s1; n2 = s2
            stripped++
        }

        print n1 rest1 > o1
        print n2 rest2 > o2
        for (i = 1; i <= 3; i++) {
            print body1[i] > o1
            print body2[i] > o2
        }
    }
}
END {
    if (!aborted)
        printf("normalizePairedReadNames: %d read pairs, %d mate-suffix stripped\n", \
               pairs, stripped) | "cat 1>&2"
}
'
