#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process checkUniqueIds {
  container = 'veupathdb/dnaseqanalysis:1.0.0'

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
  container = 'veupathdb/dnaseqanalysis:1.0.0'
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


process makeSnpFile {
  container = 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    path 'merged.vcf.gz'

  output: 
    path 'snpFile.tsv', emit: snpFile

  script:
    """
    cp merged.vcf.gz hold.vcf.gz
    gunzip hold.vcf.gz
    perl /usr/bin/makeSnpFile.pl --vcf hold.vcf --output snpFile.tsv
    """

  stub:
    """
    touch snpFile.tsv
    """
}


process processSeqVars {
  container = 'veupathdb/dnaseqanalysis:1.0.0'
  publishDir "$params.cacheFileDir", mode: "copy", pattern: "$params.cacheFile"
  publishDir "$params.outputDir", mode: "copy", pattern: 'allele.dat'
  publishDir "$params.outputDir", mode: "copy", pattern: 'product.dat'
  publishDir "$params.outputDir", mode: "copy", pattern: 'variationFeature.dat'

  input:
    path snpFile
    path cacheFile
    path undoneStrainsFile
    val  reference_strain
    path coverageDir
    path transcriptDb
    path indelDb
    path gtfFile
    path coverageComplete
    path bigwigsComplete
    path bamsComplete

  output:
    path cacheFile
    path 'variationFeature.dat', emit: variationFile
    path 'allele.dat', emit: alleleFile
    path 'product.dat', emit: productFile

  script:
    """
    set -euo pipefail

    julia /usr/bin/processSequenceVariations.jl \\
      --snp_file $snpFile \\
      --cache_file $cacheFile \\
      --undone_strains_file $undoneStrainsFile \\
      --reference_strain $reference_strain \\
      --coverage_directory $coverageDir \\
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

process addFeatureIdsToVariation {
  publishDir "$params.outputDir", mode: "copy", pattern: 'variationFeatureFinal.dat'

  input:
    path variationFile
    path gusConfig

  output:
    path 'variationFeatureFinal.dat'

  script:
    """
    set -euo pipefail
    addFeatureIdsToVariation.pl \\
         --variationFile $variationFile \\
         --gusConfig $gusConfig
    """

  stub:
    """
    touch variationFeatureFinal.dat
    """
}

process insertVariation {
  tag "plugin"

  input:
    val extDbRlsSpec
    path variationFile

  output:
    stdout

  script:
    """
    set -euo pipefail

    ga ApiCommonData::Load::Plugin::InsertVariant \\
      --extDbRlsSpec '$extDbRlsSpec' \\
      --variantFile '$variationFile' \\
      --commit

    echo "DONE"
    """

  stub:
    """
    echo "insert variation"
    """
}


process insertProduct {
  tag "plugin"

  input:
    path productFile

  output:
    stdout

  script:
    """
    set -euo pipefail

    ga ApiCommonData::Load::Plugin::InsertVariantProductSummary \\
      --variantProductFile '$productFile' \\
      --commit

    echo "DONE"
    """

  stub:
    """
    echo "insert product"
    """
}


process insertAllele {
  tag "plugin"

  input:
    path alleleFile

  output:
    stdout

  script:
    """
    set -euo pipefail

    ga ApiCommonData::Load::Plugin::InsertVariantAlleleSummary \\
      --variantAlleleFile '$alleleFile' \\
      --commit

    echo "DONE"
    """

  stub:
    """
    echo "insert allele"
    """
}


process snpEff {
  container = 'veupathdb/dnaseqanalysis:1.0.0'
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
