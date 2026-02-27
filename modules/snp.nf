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

    # Run freebayes
    freebayes \\
      -f $genomeReorderedFasta \\
      -p $params.ploidy \\
      --gvcf \\
      --gvcf-dont-use-chunk \\
      --min-coverage $params.minCoverage \\
      --min-alternate-fraction $params.freebayesMinAltFraction \\
      $resultSortedGatkBam > freebayes.vcf

    # Split into SNPs and indels for use in CNV
    bcftools view -v snps freebayes.vcf > freebayes.snps.vcf
    bcftools view -v indels freebayes.vcf > freebayes.indels.vcf

    # Compress and index split VCFs
    bgzip freebayes.snps.vcf
    tabix -fp vcf freebayes.snps.vcf.gz
    bgzip freebayes.indels.vcf
    tabix -fp vcf freebayes.indels.vcf.gz

    # Compress and index unfiltered VCF; name by sample so files are unique when merged across samples
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
    touch ${sampleName}.coverage.txt
    """
}

process makeIndelTSV {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir", pattern: "output.tsv", mode: "copy", saveAs: { filename -> "${sampleName}.indel.tsv" }

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

process bcftoolsConsensusAndMask {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple val(sampleName), path(concatVcfGz), path(concatVcfGzTbi)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)
    path(bamFile)

  output:
    tuple val(sampleName), path('cons_masked.fa')

  script:
    """
    set -euo pipefail

    bcftools consensus \\
      -I \\
      -f $genomeReorderedFasta $concatVcfGz > cons.fa

    samtools mpileup -f cons.fa -A -B ${sampleName}.bam > temp.pileup

    # Index the unmasked consensus
    samtools faidx cons.fa

    perl /usr/bin/maskGenome.pl -p temp.pileup -f cons.fa.fai -dc $params.minCoverage -o masked.fa
    fold -w 60 masked.fa > cons_masked.fa
    rm temp.pileup
    """

  stub:
    """
    touch cons_masked.fa
    """

}

process addSampleToDefline {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir", mode: "copy", saveAs: { filename -> "${sampleName}_consensus.fa.gz" }

  input:
  tuple val(sampleName), path(consFasta)

  output:
    path 'unique_ids.fa.gz'

  script:
    """
    set -euo pipefail
    addSampleToDefline.pl \\
         -i $consFasta \\
         -o unique_ids.fa \\
         -s $sampleName
    gzip unique_ids.fa
    """

  stub:
    """
    touch unique_ids.fa.gz
    """

}
