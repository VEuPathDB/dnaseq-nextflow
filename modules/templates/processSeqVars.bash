#!/usr/bin/env bash

set -euo pipefail

perl /usr/bin/processSequenceVariationsNew.pl \
  --new_sample_file snpFile.tsv \
  --cache_file cache.txt \
  --undone_strains_file undoneStrains.txt \
  --transcript_extdb_spec $transcript_extdb_spec \
  --organism_abbrev $organism_abbrev \
  --reference_strain $reference_strain  \
  --extdb_spec $extdb_spec \
  --varscan_directory varscan_directory/ \
  --gusConfigFile gusConfig.txt \
  --genome genome.fa \
  --consensus consensus.fa
