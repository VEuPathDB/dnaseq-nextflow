#!/usr/bin/env bash

set -euo pipefail

cp $consensusFasta unzipped.fa.gz;
gunzip unzipped.fa.gz;

perl /usr/bin/processSequenceVariationsNew.pl \
  --new_sample_file $snpFile \
  --cache_file $cacheFile \
  --undone_strains_file $undoneStrainsFile \
  --organism_abbrev $organism_abbrev \
  --reference_strain $reference_strain  \
  --varscan_directory $varscanDir \
  --genome $genomeFasta \
  --consensus unzipped.fa \
  --indelFile $indelFile \
  --gtfFile $gtfFile
