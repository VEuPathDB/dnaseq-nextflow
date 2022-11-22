#!/usr/bin/env bash

set -euo pipefail

cp consensus.fa.gz con.fa.gz
gunzip con.fa.gz;

if [ `grep '>' con.fa | wc -l` -eq `grep '>' con.fa | sort -u | wc -l` ];
then
     echo 'Checked all deflines... no repeats';
else
     echo 'ERROR:  Defines in .fa files must be unique';  exit 125;
fi
