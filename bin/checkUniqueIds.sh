#!/usr/bin/env bash

set -euo pipefail

mkdir FASTAS;
cp $1/*.fa.gz FASTAS;
gunzip FASTAS/*.gz;

if [ `grep '>' FASTAS/*.fa | wc -l` -eq `grep '>' FASTAS/*.fa | sort -u | wc -l` ];
then
     echo 'Checked all deflines... no repeats';
else
     echo 'ERROR:  Defines in .fa files must be unique';  exit 125;
fi
