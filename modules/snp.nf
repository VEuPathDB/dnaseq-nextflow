#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process freebayes {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple val(sampleName), path(resultSortedGatkBam), path(resultSortedGatkBamIndex)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    tuple val(sampleName), path("${sampleName}.vcf.gz"), path("${sampleName}.vcf.gz.tbi"), path('freebayes.snps.vcf.gz'), path('freebayes.snps.vcf.gz.tbi'), path('freebayes.indels.vcf.gz'), path('freebayes.indels.vcf.gz.tbi'), emit: vcf_files

  script:
    """
    set -euo pipefail

    minAltFraction=\$([ "$params.ploidy" -eq 1 ] && echo "0.8" || echo "0.3")
    freebayes \\
      -f $genomeReorderedFasta \\
      -p $params.ploidy \\
      --min-coverage $params.minCoverage \\
      --min-alternate-fraction \$minAltFraction \\
      $resultSortedGatkBam | vcfallelicprimitives -kg | bcftools sort > freebayes.vcf

    bcftools view -v snps freebayes.vcf > freebayes.snps.vcf
    bcftools norm -m- freebayes.vcf | \\
      bcftools view --include 'strlen(ALT)!=strlen(REF) && ALT!~"^<"' > freebayes.indels.vcf

    bgzip freebayes.snps.vcf
    tabix -fp vcf freebayes.snps.vcf.gz
    bgzip freebayes.indels.vcf
    tabix -fp vcf freebayes.indels.vcf.gz

    bgzip freebayes.vcf
    mv freebayes.vcf.gz ${sampleName}.vcf.gz
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


process mergeVcfs {
  container 'biocontainers/bcftools:v1.9-1-deb_cv1'

  input:
    val vcfCount
    tuple path(samplevcfzip), path(samplevcfzipindex), val(key)

  output:
    path('result.vcf.gz')

  script:
    """
    set -euo pipefail

    if [ $vcfCount -gt 1 ]; then
        bcftools merge -o result.vcf.gz -O z *.vcf.gz
    else
        cp *.vcf.gz result.vcf.gz
    fi
    """

  stub:
    """
    touch result.vcf.gz
    """

}

process makeMergedVariantIndex {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir", mode: "copy"

  input:
    path(resultVcfGz)

  output:
    tuple path('result.vcf.gz'), path('result.vcf.gz.tbi')

  script:
    """
    set -euo pipefail
    cp $resultVcfGz hold.vcf.gz
    gunzip hold.vcf.gz
    sed -i 's/\\%//g' hold.vcf
    bgzip hold.vcf
    mv hold.vcf.gz result.vcf.gz
    tabix -fp vcf result.vcf.gz
    """

  stub:
    """
    touch result.vcf.gz
    touch result.vcf.gz.tbi
    """

}

process bcftoolsMpileupGvcf {
  container 'biocontainers/bcftools:v1.9-1-deb_cv1'

  input:
    tuple val(sampleName), path(bamFile), path(bamIndex)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    tuple val(sampleName), path("${sampleName}.g.vcf.gz"), path("${sampleName}.g.vcf.gz.tbi")

  script:
    """
    set -euo pipefail
    bcftools mpileup \\
      --gvcf $params.minCoverage \\
      -f $genomeReorderedFasta \\
      $bamFile \\
    | bcftools call \\
      -m \\
      --gvcf $params.minCoverage \\
      -O z \\
      -o ${sampleName}.g.vcf.gz
    bcftools index -t ${sampleName}.g.vcf.gz
    """

  stub:
    """
    touch ${sampleName}.g.vcf.gz
    touch ${sampleName}.g.vcf.gz.tbi
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

  publishDir "$params.outputDir", mode: "copy", saveAs: { "${sampleName}_consensus.fa.gz" }

  input:
    tuple val(sampleName), path(vcfGz), path(vcfGzTbi), path(coverageBed)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    path 'consensus.fa.gz'

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
    """

  stub:
    """
    touch consensus.fa.gz
    """
}
