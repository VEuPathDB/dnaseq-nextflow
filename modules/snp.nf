#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process runFreebayes {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple val(sampleName), path(resultSortedGatkBam), path(resultSortedGatkBamIndex)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    tuple val(sampleName), path("${sampleName}.vcf.gz"), path("${sampleName}.vcf.gz.tbi")

  script:
    """
    set -euo pipefail

    minAltFraction=\$([ "$params.ploidy" -eq 1 ] && echo "0.8" || echo "0.3")
    freebayes \\
      -f $genomeReorderedFasta \\
      -p $params.ploidy \\
      --min-coverage $params.minCoverage \\
      --min-alternate-fraction \$minAltFraction \\
      $resultSortedGatkBam | bcftools norm -f $genomeReorderedFasta -a | mergeVariantsByLocation.py | bcftools sort > freebayes.vcf

    bgzip freebayes.vcf
    mv freebayes.vcf.gz ${sampleName}.vcf.gz
    tabix -fp vcf ${sampleName}.vcf.gz
    """

  stub:
    """
    touch ${sampleName}.vcf.gz
    touch ${sampleName}.vcf.gz.tbi
    """
}

process filterAndSplitVcf {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple val(sampleName), path(rawVcfGz), path(rawVcfGzTbi)

  output:
    tuple val(sampleName), path("${sampleName}.vcf.gz"), path("${sampleName}.vcf.gz.tbi"), path('freebayes.snps.vcf.gz'), path('freebayes.snps.vcf.gz.tbi'), path('freebayes.indels.vcf.gz'), path('freebayes.indels.vcf.gz.tbi'), emit: vcf_files

  script:
    """
    set -euo pipefail

    bcftools view -e 'GT="0/0"' $rawVcfGz | bcftools filter -e 'RPP > 20' > filtered.vcf

    bcftools view -v snps filtered.vcf > freebayes.snps.vcf
    bcftools norm -m- filtered.vcf | \\
      bcftools view --include 'strlen(ALT)!=strlen(REF) && ALT!~"^<" && ALT!="*"' > freebayes.indels.vcf

    bgzip freebayes.snps.vcf
    tabix -fp vcf freebayes.snps.vcf.gz
    bgzip freebayes.indels.vcf
    tabix -fp vcf freebayes.indels.vcf.gz

    bgzip filtered.vcf
    mv filtered.vcf.gz ${sampleName}.vcf.gz
    tabix -fp vcf ${sampleName}.vcf.gz
    """

  stub:
    """
    touch ${sampleName}.vcf.gz
    touch ${sampleName}.vcf.gz.tbi
    touch freebayes.snps.vcf.gz
    touch freebayes.snps.vcf.gz.tbi
    touch freebayes.indels.vcf.gz
    touch freebayes.indels.vcf.gz.tbi
    """
}

process makeIndelTSV {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple val(sampleName), path(indelsVcfGz)

  output:
    path('output.tsv')

  script:
    """
    set -euo pipefail
    findValues.pl \\
       -i $indelsVcfGz \\
       -s ${sampleName} \\
       -o output.tsv
    """

  stub:
    """
    touch output.tsv
    """

}

process sanitizeVcf {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir", mode: "copy"

  input:
    tuple val(sampleName), path(vcfGz), path(vcfGzTbi)

  output:
    tuple val(sampleName), path("${sampleName}.vcf.gz"), path("${sampleName}.vcf.gz.tbi")

  script:
    """
    set -euo pipefail
    gunzip -c $vcfGz | sed 's/%//g' | bgzip > ${sampleName}.vcf.gz
    tabix -fp vcf ${sampleName}.vcf.gz
    """

  stub:
    """
    touch ${sampleName}.vcf.gz
    touch ${sampleName}.vcf.gz.tbi
    """
}

process makeCoverageBed {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir", mode: "copy"

  input:
    tuple val(sampleName), path(bamFile), path(bamIndex)

  output:
    tuple val(sampleName), path("${sampleName}.coverage.bed.gz")

  script:
    """
    set -euo pipefail
    bedtools genomecov -ibam $bamFile -bga \\
      | awk -v mc=$params.minCoverage '\$4 >= mc' \\
      | bedtools merge -c 4 -o mean \\
      | bgzip > ${sampleName}.coverage.bed.gz
    """

  stub:
    """
    touch ${sampleName}.coverage.bed.gz
    """
}


process makeConsensusFromVcfAndBed {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir", mode: "copy"

  input:
    tuple val(sampleName), path(vcfGz), path(vcfGzTbi), path(coverageBed)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    path "${sampleName}_consensus.fa.gz"

  script:
    """
    set -euo pipefail
    makeConsensusFastaFromVcfAndBed.py \\
      --vcf $vcfGz \\
      --bed $coverageBed \\
      --ref $genomeReorderedFasta \\
      --fai $genomeReorderedFastaIndex \\
      --output consensus.fa
    bgzip consensus.fa
    mv consensus.fa.gz ${sampleName}_consensus.fa.gz
    """

  stub:
    """
    touch consensus.fa.gz
    """
}
