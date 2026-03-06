#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process checkUniqueIds {
  container 'veupathdb/dnaseqanalysis:1.0.0'
  input:
    path 'consensus.fa.gz'

  output:
    stdout

  script:
    """
    set -euo pipefail
    checkUniqueIds.sh
    """

  stub:
    """
    touch stdout
    """
}


process mergeVcfs {
  container 'veupathdb/dnaseqanalysis:1.0.0'
  publishDir "$params.outputDir", mode: "copy", pattern: 'merged.vcf.gz'

  input:
    path '*.vcf.gz'

  output:
    path 'merged.vcf.gz'

  script:
    """
    set -euo pipefail

    for i in *.vcf.gz; do cp \$i \$i.tmp.vcf.gz; gunzip \$i.tmp.vcf.gz; bgzip \$i.tmp.vcf; cp \$i.tmp.vcf.gz \$i; rm \$i.tmp.vcf.gz; tabix \$i; done
    bcftools merge \\
          -o merged.vcf.gz \\
          -O z *.vcf.gz
    cp merged.vcf.gz merge.vcf.gz
    gunzip merge.vcf.gz
    sed -i 's/\\%//g' merge.vcf
    mv merge.vcf merged.vcf
    """

  stub:
    """
    touch merged.vcf.gz
    touch merged.vcf
    """
}


process makeCodingData {
  container 'veupathdb/shortreadaligner:1.0.0'

  input:
    path fastas
    path genomicIndelDb
    path gtfFile
    path genomeFastaFile

  output:
    path 'codingSequences.db', emit: codingSequencesDb
    path 'codingIndels.db',    emit: codingIndelsDb

  script:
    """
    set -euo pipefail
    makeCodingData.jl \\
      --genomic_indel_db $genomicIndelDb \\
      --gtf_file $gtfFile \\
      --genome_fasta $genomeFastaFile \\
      --cds_db_out codingSequences.db \\
      --indels_db_out codingIndels.db
    """

  stub:
    """
    touch codingSequences.db
    touch codingIndels.db
    """
}


process makeGenomicIndelDb {
  container 'veupathdb/shortreadaligner:1.0.0'

  input:
    path 'indels.tsv'

  output:
    path 'genomicIndels.db'

  script:
    """
    set -euo pipefail
    sqlite3 genomicIndels.db <<'SQL'
CREATE TABLE genomic_indels (
  strain      TEXT    NOT NULL,
  sequence_id TEXT    NOT NULL,
  position    INTEGER NOT NULL,
  shift       INTEGER NOT NULL
);
.separator "\\t"
.import indels.tsv genomic_indels
CREATE INDEX idx_genomic_indels ON genomic_indels(sequence_id, strain, position);
SQL
    """

  stub:
    """
    touch genomicIndels.db
    """
}


process processSeqVars {
  container 'veupathdb/dnaseqanalysis:1.0.0'
  publishDir "$params.outputDir", mode: "copy", pattern: "$params.vcfCacheFile"
  publishDir "$params.outputDir", mode: "copy", pattern: 'allele.dat'
  publishDir "$params.outputDir", mode: "copy", pattern: 'product.dat'
  publishDir "$params.outputDir", mode: "copy", pattern: 'variationFeature.dat'

  input:
    path vcfFile
    path cacheFile
    path undoneStrainsFile
    val  reference_strain
    path transcriptDb
    path indelDb
    path gtfFile

  output:
    path cacheFile
    path 'variationFeature.dat', emit: variationFile
    path 'allele.dat', emit: alleleFile
    path 'product.dat', emit: productFile

  script:
    """
    set -euo pipefail

    julia /usr/bin/processSequenceVariations.jl \\
      --vcf_file $vcfFile \\
      --cache_file $cacheFile \\
      --undone_strains_file $undoneStrainsFile \\
      --reference_strain $reference_strain \\
      --transcript_db $transcriptDb \\
      --indel_db $indelDb \\
      --gtf_file $gtfFile

    mv snpFeature.dat variationFeature.dat
    """

  stub:
    """
    touch cache.txt
    touch snpFeature.dat
    touch allele.dat
    touch product.dat
    """
}


process snpEff {
  container 'veupathdb/dnaseqanalysis:1.0.0'
  publishDir "$params.outputDir", mode: "copy"

  input:
    path mergedVcf
    path genesGtf
    path sequencesFa

  output:
    path 'merged.ann.vcf'

  script:
    """
    set -euo pipefail
    mkdir genome
    mv $genesGtf genome/genes.gtf
    mv $sequencesFa genome/sequences.fa
    gzip -f genome/sequences.fa
    cp /usr/bin/snpEff/snpEff.config .
    java -jar /usr/bin/snpEff/snpEff.jar build -gtf22 -noCheckCds -noCheckProtein -v genome
    java -Xmx4g -jar /usr/bin/snpEff/snpEff.jar genome $mergedVcf > merged.ann.vcf
    """

  stub:
    """
    touch merged.ann.vcf
    """
}
