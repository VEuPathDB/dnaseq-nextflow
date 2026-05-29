#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process checkUniqueIds {
  container 'veupathdb/dnaseqanalysis:1.0.1'
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
  container 'veupathdb/dnaseqanalysis:1.0.1'

  input:
    path "*.vcf.gz"

  output:
    path 'merged.vcf.gz'

  script:
    """
    set -euo pipefail

    for vcf in *.vcf.gz; do bcftools index --tbi \$vcf; done
    bcftools merge --merge all -O z -o merged.vcf.gz *.vcf.gz
    """

  stub:
    """
    touch merged.vcf.gz
    """
}



process mergeCoverageBeds {
  container 'veupathdb/dnaseqanalysis:1.0.1'

  input:
    path coverageBeds

  output:
    path 'coverage.tsv'

  script:
    """
    set -euo pipefail
    files=( $coverageBeds )
    names=()
    for f in "\${files[@]}"; do
      names+=( "\$(basename "\$f" .coverage.bed.gz)" )
    done
    header="chrom\tstart\tend\t\$(IFS='\t'; echo "\${names[*]}")"
    echo -e "\$header" > coverage.tsv
    bedtools unionbedg -names "\${names[@]}" -filler 0 -i "\${files[@]}" >> coverage.tsv
    """

  stub:
    """
    touch coverage.tsv
    """
}


process makeCodingData {
  container 'veupathdb/dnaseqanalysis:1.0.1'

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
  container 'veupathdb/dnaseqanalysis:1.0.1'

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
  container 'veupathdb/dnaseqanalysis:1.0.1'

  publishDir "$params.outputDir", mode: "copy", pattern: 'allele.dat'
  publishDir "$params.outputDir", mode: "copy", pattern: 'transcript_product.dat'
  publishDir "$params.outputDir", mode: "copy", pattern: 'variationFeature.dat'
  publishDir "$params.outputDir", mode: "copy", pattern: 'sample.dat'
  publishDir "$params.outputDir", mode: "copy", pattern: 'hsss_readFreq20'
  publishDir "$params.outputDir", mode: "copy", pattern: 'hsss_readFreq40'
  publishDir "$params.outputDir", mode: "copy", pattern: 'hsss_readFreq60'
  publishDir "$params.outputDir", mode: "copy", pattern: 'hsss_readFreq80'

  input:
    path vcfFile
    path cacheFile  // previous run's transcript_product.dat (was cache.tsv)
    path undoneStrainsFile
    val  reference_strain
    path transcriptDb
    path indelDb
    path gtfFile
    path coverageFile

  output:
    path 'output.vcf.gz', emit: outputVcf
    path 'output.vcf.gz.tbi', emit: outputVcfIndex
    path 'transcript_product.dat', emit: transcriptProductFile
    path 'variationFeature.dat', emit: variationFile
    path 'allele.dat', emit: alleleFile
    path 'sample.dat', emit: sampleFile
    path 'hsss_readFreq20', emit: hsssReadFreq20
    path 'hsss_readFreq40', emit: hsssReadFreq40
    path 'hsss_readFreq60', emit: hsssReadFreq60
    path 'hsss_readFreq80', emit: hsssReadFreq80

  script:
    """
    set -euo pipefail

    processSequenceVariations.jl \\
      --vcf_file $vcfFile \\
      --cache_file $cacheFile \\
      --undone_strains_file $undoneStrainsFile \\
      --reference_strain $reference_strain \\
      --transcript_db $transcriptDb \\
      --indel_db $indelDb \\
      --gtf_file $gtfFile \\
      --coverage_file $coverageFile \\
      --ploidy $params.ploidy \\
      --output_vcf output.vcf

    mv snpFeature.dat variationFeature.dat
    bgzip output.vcf
    tabix -p vcf output.vcf.gz
    """

  stub:
    """
    touch output.vcf.gz
    touch output.vcf.gz.tbi
    touch transcript_product.dat
    touch variationFeature.dat
    touch allele.dat
    touch sample.dat
    mkdir hsss_readFreq20
    mkdir hsss_readFreq40
    mkdir hsss_readFreq60
    mkdir hsss_readFreq80
    """
}


process snpEff {
  container 'veupathdb/dnaseqanalysis:1.0.1'
  publishDir "$params.outputDir", mode: "copy"

  input:
    path mergedVcf
    path genesGtf
    path sequencesFa

  output:
    path 'merged.ann.vcf.gz',     emit: annotatedVcf
    path 'merged.ann.vcf.gz.tbi', emit: annotatedVcfIndex

  script:
    """
    set -euo pipefail
    mkdir genome
    mv $genesGtf genome/genes.gtf
    mv $sequencesFa genome/sequences.fa
    gzip -f genome/sequences.fa
    cp /usr/bin/snpEff/snpEff.config .
    java -jar /usr/bin/snpEff/snpEff.jar build -gtf22 -noCheckCds -noCheckProtein -v genome
    java -Xmx4g -jar /usr/bin/snpEff/snpEff.jar -no-downstream -no-upstream genome $mergedVcf > merged.ann.vcf
    bgzip merged.ann.vcf
    tabix -p vcf merged.ann.vcf.gz
    """

  stub:
    """
    touch merged.ann.vcf.gz merged.ann.vcf.gz.tbi
    """
}


process parseSnpEffAnnotations {
  container 'veupathdb/shortreadaligner:1.0.0'
  publishDir "$params.outputDir", mode: "copy", pattern: 'snpeff.dat'

  input:
    path annVcf

  output:
    path 'snpeff.dat', emit: snpeffFile

  script:
    """
    set -euo pipefail
    parseSnpEffAnnotations.py --vcf $annVcf --output snpeff.dat
    """

  stub:
    """
    touch snpeff.dat
    """
}
