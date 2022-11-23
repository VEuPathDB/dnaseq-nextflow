#!/usr/bin/env bash

set -euo pipefail

cp consensus.fa.gz unzipped.fa.gz;
cp genome.fa.gz hold.fa.gz

gunzip unzipped.fa.gz;
gunzip hold.fa.gz;

processSequenceVariationsNew.pl \
  --new_sample_file snpFile.tsv \
  --cache_file cache.txt \
  --undone_strains_file undoneStrains.txt \
  --transcript_extdb_spec $genomeExtDbRlsSpec \
  --organism_abbrev $organism_abbrev \
  --reference_strain $reference_strain  \
  --extdb_spec $extdb_spec \
  --varscan_directory varscan_directory/ \
  --gusConfigFile gusConfig.txt \
  --genome hold.fa \
  --consensus unzipped.fa
