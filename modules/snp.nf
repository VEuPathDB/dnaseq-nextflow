#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process freebayes {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple val(sampleName), path(resultSortedGatkBam), path(resultSortedGatkBamIndex)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    tuple val(sampleName), path("${sampleName}.vcf.gz"), path("${sampleName}.vcf.gz.tbi"), path('freebayes.snps.vcf.gz'), path('freebayes.snps.vcf.gz.tbi'), path('freebayes.indels.vcf.gz'), path('freebayes.indels.vcf.gz.tbi'), path("${sampleName}.g.vcf.gz"), path("${sampleName}.g.vcf.gz.tbi"), emit: vcf_files

  script:
    """
    set -euo pipefail

    # Run freebayes with --gvcf to include reference blocks
    minAltFraction=\$([ "$params.ploidy" -eq 1 ] && echo "0.8" || echo "0.3")
    freebayes \\
      -f $genomeReorderedFasta \\
      -p $params.ploidy \\
      --min-coverage $params.minCoverage \\
      --min-alternate-fraction \$minAltFraction \\
      --gvcf \\
      $resultSortedGatkBam | bcftools sort > freebayes.g.vcf

    # Extract variant sites only (exclude reference blocks where ALT=<*>)
    bcftools view -e 'ALT[0]="<*>"' freebayes.g.vcf > freebayes.vcf

    # Split into SNPs and indels for use in CNV
    bcftools view -v snps freebayes.vcf > freebayes.snps.vcf
    # Split multi-allelic records, then select indels by length difference (catches complex variants
    # missed by -v indels) while excluding symbolic alleles like <*>
    bcftools norm -m- freebayes.vcf | \
      bcftools view --include 'strlen(ALT)!=strlen(REF) && ALT!~"^<"' > freebayes.indels.vcf

    # Compress and index split VCFs
    bgzip freebayes.snps.vcf
    tabix -fp vcf freebayes.snps.vcf.gz
    bgzip freebayes.indels.vcf
    tabix -fp vcf freebayes.indels.vcf.gz

    # Compress and index unfiltered VCF; name by sample so files are unique when merged across samples
    bgzip freebayes.vcf
    mv freebayes.vcf.gz ${sampleName}.vcf.gz
    tabix -fp vcf ${sampleName}.vcf.gz

    # Compress and index GVCF
    bgzip freebayes.g.vcf
    mv freebayes.g.vcf.gz ${sampleName}.g.vcf.gz
    tabix -fp vcf ${sampleName}.g.vcf.gz
    """

  stub:
    """
    touch ${sampleName}.vcf.gz
    touch ${sampleName}.vcf.gz.tbi
    touch freebayes.snps.vcf.gz
    touch freebayes.snps.vcf.gz.tbi
    touch freebayes.indels.vcf.gz
    touch freebayes.indels.vcf.gz.tbi
    touch ${sampleName}.g.vcf.gz
    touch ${sampleName}.g.vcf.gz.tbi
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

process mergeGvcfs {
  container 'biocontainers/bcftools:v1.9-1-deb_cv1'

  publishDir "$params.outputDir", mode: "copy"

  input:
    val gvcfCount
    tuple path(gvcfFiles), path(gvcfIndexes), val(key)

  output:
    tuple path('coverage.g.vcf.gz'), path('coverage.g.vcf.gz.tbi')

  script:
    """
    set -euo pipefail
    if [ $gvcfCount -gt 1 ]; then
    bcftools merge \
        --merge all \
        --output-type z \
        --output coverage.g.vcf.gz \
        *.g.vcf.gz
    else
        cp *.g.vcf.gz coverage.g.vcf.gz
    fi
    bcftools index -t coverage.g.vcf.gz
    """

  stub:
    """
    touch coverage.g.vcf.gz
    touch coverage.g.vcf.gz.tbi
    """
}

process splitGvcfAtZeroCoverage {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple val(regionKey), path(gvcfGz, stageAs: 'input.g.vcf.gz'), path(gvcfGzTbi, stageAs: 'input.g.vcf.gz.tbi'), path(zeroCovBed)
    tuple path(genomeFasta), path(genomeFastaIndex)

  output:
    tuple val(regionKey), path("${regionKey}.split.g.vcf.gz"), path("${regionKey}.split.g.vcf.gz.tbi")

  script:
    """
    set -euo pipefail
    splitGvcfAtZeroCoverage.py \\
      --gvcf input.g.vcf.gz \\
      --zero-cov-bed $zeroCovBed \\
      --ref $genomeFasta \\
      --output /dev/stdout \\
    | bcftools sort -O z -o ${regionKey}.split.g.vcf.gz
    bcftools index -t ${regionKey}.split.g.vcf.gz
    """

  stub:
    """
    touch ${regionKey}.split.g.vcf.gz ${regionKey}.split.g.vcf.gz.tbi
    """
}

process freebayesMultiSample {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    path bamFiles
    path bamIndexes
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)
    val regionLine

  output:
    tuple val(regionKey), path("${regionKey}.g.vcf.gz"), path("${regionKey}.g.vcf.gz.tbi"), emit: gvcf

  script:
    def fields    = regionLine.tokenize('\t')
    def chrom     = fields[0]
    def start     = fields[1].toLong()
    def end       = fields[2].toLong()
    regionKey     = "${chrom}_${start}_${end}"
    def regionStr = "${chrom}:${start}-${end}"
    """
    set -euo pipefail
    ls *.bam > bam_list.txt
    minAltFraction=\$([ "$params.ploidy" -eq 1 ] && echo "0.8" || echo "0.3")
    freebayes \\
      -f $genomeReorderedFasta \\
      -p $params.ploidy \\
      --min-coverage $params.minCoverage \\
      --min-alternate-fraction \$minAltFraction \\
      --gvcf \\
      --region ${regionStr} \\
      --bam-list bam_list.txt \\
    | bcftools sort -O z -o ${regionKey}.g.vcf.gz
    bcftools index -t ${regionKey}.g.vcf.gz
    """

  stub:
    def fields = regionLine.tokenize('\t')
    def chrom  = fields[0]
    def start  = fields[1].toLong()
    def end    = fields[2].toLong()
    regionKey  = "${chrom}_${start}_${end}"
    """
    touch ${regionKey}.g.vcf.gz ${regionKey}.g.vcf.gz.tbi
    """
}

process makeRegionBed {
  container 'veupathdb/shortreadaligner:1.0.0'

  input:
    tuple path(genomeFasta), path(genomeFastaIndex)

  output:
    path 'regions.bed'

  script:
    """
    set -euo pipefail
    bedtools makewindows -g $genomeFastaIndex -w $params.chunkSize > regions.bed
    """

  stub:
    """
    printf 'chr1\\t0\\t1000000\\n' > regions.bed
    printf 'chr1\\t1000000\\t2000000\\n' >> regions.bed
    """
}

process makeMultiSampleZeroCoverageBed {
  container 'veupathdb/shortreadaligner:1.0.0'

  input:
    path bamFiles
    path bamIndexes
    val regionLine

  output:
    tuple val(regionKey), path('zero_cov.bed')

  script:
    def fields = regionLine.tokenize('\t')
    def chrom  = fields[0]
    def start  = fields[1].toLong()
    def end    = fields[2].toLong()
    def startPlus1 = start + 1
    regionKey  = "${chrom}_${start}_${end}"
    """
    set -euo pipefail
    for bam in *.bam; do
      samtools view -b -h \$bam ${chrom}:${startPlus1}-${end} > region_\${bam%.bam}.bam
      samtools index region_\${bam%.bam}.bam
      bedtools genomecov -ibam region_\${bam%.bam}.bam -bga | \\
        awk '\$4 == 0 {print \$1 "\\t" \$2 "\\t" \$3}' >> all_zero.bed
    done
    if [ -s all_zero.bed ]; then
      bedtools sort -i all_zero.bed | bedtools merge > zero_cov.bed
    else
      touch zero_cov.bed
    fi
    """

  stub:
    def fields = regionLine.tokenize('\t')
    def chrom  = fields[0]
    def start  = fields[1].toLong()
    def end    = fields[2].toLong()
    regionKey  = "${chrom}_${start}_${end}"
    """
    touch zero_cov.bed
    """
}

process makeConsensusFromGvcf {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir", mode: "copy", saveAs: { "${sampleName}_consensus.fa.gz" }

  input:
    tuple val(sampleName), path(gvcfGz), path(gvcfGzTbi)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    path 'consensus.fa.gz'

  script:
    """
    set -euo pipefail
    makeConsensusFastaFromGvcf.py \\
      --gvcf $gvcfGz \\
      --ref $genomeReorderedFasta \\
      --fai $genomeReorderedFastaIndex \\
      --min-coverage $params.minCoverage \\
      --output consensus.fa
    bgzip consensus.fa
    """

  stub:
    """
    touch consensus.fa.gz
    """

}

