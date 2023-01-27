#!/usr/bin/env bash

set -euo pipefail

cp $consensusFasta unzipped.fa.gz;
gunzip unzipped.fa.gz;

cp $snpFile snp.tsv
cp $cacheFile cf.txt
cp $undoneStrainsFile us.txt
mkdir varCons
cp $varscanDir/* varCons
cp $genomeFasta genome.fasta
cp $indelFile if.tsv
cp $gtfFile gtf.gtf

perl /usr/bin/processSequenceVariationsNew.pl \
  --new_sample_file snp.tsv \
  --cache_file cf.txt \
  --undone_strains_file us.txt \
  --organism_abbrev $organism_abbrev \
  --reference_strain $reference_strain  \
  --varscan_directory ./varCons \
  --genome genome.fasta \
  --consensus unzipped.fa \
  --indelFile if.tsv \
  --gtfFile gtf.gtf
