#!/usr/bin/env bash

set -euo pipefail

FASTA_FILES=$fastaDir/*.fa

if [ "$(grep '>' $FASTA_FILES | wc -l )" -eq "$(grep '>' $FASTA_FILES | sort -u | wc -l )" ];
then
     echo 'Checked all deflines... no repeats';
else
     echo 'ERROR:  Defines in .fa files must be unique';  exit 125;
fi
