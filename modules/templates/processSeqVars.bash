#!/usr/bin/env bash

set -euo pipefail

cp consensus.fa.gz unzipped.fa.gz;
gunzip unzipped.fa.gz;

processSequenceVariationsNew.pl \
  --new_sample_file snpFile.tsv \
  --cache_file cache.txt \
  --undone_strains_file undoneStrains.txt \
  --transcript_extdb_spec "$genomeExtDbRlsSpec" \
  --organism_abbrev $organism_abbrev \
  --reference_strain $reference_strain  \
  --extdb_spec "$extdb_spec" \
  --varscan_directory varscan_directory/ \
  --gusConfigFile gusConfig.txt \
  --genome genome.fa \
  --consensus unzipped.fa
