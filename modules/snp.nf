#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process runFreebayes {
  container 'veupathdb/dnaseqanalysis:1.0.1'

  input:
    tuple val(sampleName), path(resultSortedGatkBam), path(resultSortedGatkBamIndex)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)
    path genomeRepeatsBed

  output:
    tuple val(sampleName), path("${sampleName}.vcf.gz"), path("${sampleName}.vcf.gz.tbi")

  script:
    def targetsArg = genomeRepeatsBed.name != 'NO_FILE' ? "--targets targets.bed" : ""
    """
    set -euo pipefail

    if [ -n "$targetsArg" ]; then
      awk 'OFS="\\t"{print \$1,\$2}' ${genomeReorderedFasta}.fai > genome.sizes
      bedtools complement -i $genomeRepeatsBed -g genome.sizes > targets.bed
    fi

    minAltFraction=\$([ "$params.ploidy" -eq 1 ] && echo "0.8" || echo "0.3")
    freebayes \\
      -f $genomeReorderedFasta \\
      -p $params.ploidy \\
      --min-coverage $params.minCoverage \\
      --min-alternate-fraction \$minAltFraction \\
      $targetsArg \\
      $resultSortedGatkBam | bcftools norm -f $genomeReorderedFasta | bcftools sort > freebayes.vcf

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

// filterAndSplitVcf takes the raw complex VCF from runFreebayes and:
//   - sanitizes (strips FreeBayes % artifacts)
//   - produces consensus_input.vcf.gz: complex filtered, used for consensus FASTA generation
//   - produces freebayes.indels.vcf.gz: indels extracted from complex VCF
//   - produces ${sampleName}.vcf.gz: atomized+filtered, published for mergeExperiments
//   - produces freebayes.snps.vcf.gz: SNPs extracted from atomized VCF
process filterAndSplitVcf {
  container 'veupathdb/dnaseqanalysis:1.0.1'

  publishDir "$params.outputDir/${sampleName}", mode: "copy", pattern: "${sampleName}.vcf.gz*"

  input:
    tuple val(sampleName), path(rawVcfGz), path(rawVcfGzTbi)

  output:
    tuple val(sampleName),
      path("${sampleName}.vcf.gz"), path("${sampleName}.vcf.gz.tbi"),
      path('freebayes.snps.vcf.gz'), path('freebayes.snps.vcf.gz.tbi'),
      path('freebayes.indels.vcf.gz'), path('freebayes.indels.vcf.gz.tbi'),
      path('consensus_input.vcf.gz'), path('consensus_input.vcf.gz.tbi'),
      emit: vcf_files

  script:
    """
    set -euo pipefail

    # Sanitize: strip FreeBayes % artifacts
    bgzip -d -c $rawVcfGz | sed 's/%//g' > sanitized.vcf

    # Filter: drop hom-ref calls and high strand-bias
    bcftools view -e 'GT="0/0"' sanitized.vcf | bcftools filter -e 'RPP > 20' > complex_filtered.vcf

    # Consensus input: complex filtered (preserves complex variants for accurate FASTA)
    bgzip -c complex_filtered.vcf > consensus_input.vcf.gz
    tabix -fp vcf consensus_input.vcf.gz

    # Indels from complex VCF (correctly represents deletion component of complex variants)
    bcftools norm -m- complex_filtered.vcf | \\
      bcftools view --include 'strlen(ALT)!=strlen(REF) && ALT!~"^<" && ALT!="*"' > freebayes.indels.vcf
    bgzip freebayes.indels.vcf
    tabix -fp vcf freebayes.indels.vcf.gz

    # Atomize for mergeExperiments (single line per SNP across samples)
    bcftools norm -a complex_filtered.vcf | mergeVariantsByLocation.py | bcftools sort > atomized_filtered.vcf

    # SNPs from atomized
    bcftools view -v snps atomized_filtered.vcf > freebayes.snps.vcf
    bgzip freebayes.snps.vcf
    tabix -fp vcf freebayes.snps.vcf.gz

    # Published per-sample VCF: atomized, for mergeExperiments
    bgzip atomized_filtered.vcf
    mv atomized_filtered.vcf.gz ${sampleName}.vcf.gz
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
    touch consensus_input.vcf.gz
    touch consensus_input.vcf.gz.tbi
    """
}

process makeIndelTSV {
  container 'veupathdb/dnaseqanalysis:1.0.1'

  publishDir "${params.outputDir}/${sampleName}", mode: "copy"

  input:
    tuple val(sampleName), path(indelsVcfGz)

  output:
    path("${sampleName}_indels.tsv")

  script:
    """
    set -euo pipefail
    findValues.pl \\
       -i $indelsVcfGz \\
       -s ${sampleName} \\
       -o ${sampleName}_indels.tsv
    """

  stub:
    """
    touch ${sampleName}_indels.tsv
    """

}


process makeCoverageBed {
  container 'veupathdb/dnaseqanalysis:1.0.1'

  publishDir "$params.outputDir/${sampleName}", mode: "copy"

  input:
    tuple val(sampleName), path(bamFile), path(bamIndex)

  output:
    tuple val(sampleName), path("${sampleName}_coverage.bed.gz")

  script:
    """
    set -euo pipefail
    bedtools genomecov -ibam $bamFile -bga \\
      | awk -v mc=$params.minCoverage '\$4 >= mc' \\
      | bedtools merge -c 4 -o mean \\
      | bgzip > ${sampleName}_coverage.bed.gz
    """

  stub:
    """
    touch ${sampleName}_coverage.bed.gz
    """
}


process makeConsensusFromVcfAndBed {
  container 'veupathdb/dnaseqanalysis:1.0.1'

  publishDir "$params.outputDir/${sampleName}", mode: "copy"

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
    touch ${sampleName}_consensus.fa.gz
    """
}
