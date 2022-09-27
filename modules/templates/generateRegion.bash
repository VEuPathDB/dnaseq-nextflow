#!/usr/bin/env bash

set -euo pipefail

grep ">" input.fasta > temp.fa
DEFLINE=\$(sed 's/>//' temp.fa)
/usr/bin/perl makepositionarraycoding.pl \
                --test_file shifted.txt  \
                --sequence_id \$DEFLINE \
		--region_file region.txt
