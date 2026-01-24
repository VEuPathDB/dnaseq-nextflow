#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process mpileup {
  container = 'veupathdb/shortreadaligner:1.0.0'

  publishDir "$params.outputDir", pattern: "result.pileup", mode: "copy", saveAs: { filename -> "${sampleName}.result.pileup" }

  input:
    tuple val(sampleName), path (resultSortedGatkBam), path(resultSortedGatkBamIndex)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    tuple val(sampleName), path('result.pileup')

  script:
    """
    set -euo pipefail
    samtools mpileup \\
      -A \\
      -f $genomeReorderedFasta \\
      -B $resultSortedGatkBam > result.pileup 2>pileup.err
    """

  stub:
    """
    touch result.pileup
    """

}

process varscan {
  container = 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir/varscanCons", pattern: "${sampleName}.coverage.txt", mode: "copy"

  input:
    tuple val(sampleName), path (resultSortedGatkBam), path(resultSortedGatkBamIndex), path(resultPileup)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi'), path('varscan.indels.vcf.gz'), path('varscan.indels.vcf.gz.tbi'), path('genome_masked.fa'), emit: vcf_files
    path "${sampleName}.coverage.txt"

  script:
    """
    set -euo pipefail

    JARPATH="/usr/local/VarScan.jar"
    echo $sampleName >vcf_sample_name

    java -jar \$JARPATH mpileup2snp $resultPileup --vcf-sample-list vcf_sample_name --output-vcf 1 --p-value $params.varscanPValue --min-coverage $params.minCoverage --min-var-freq $params.varscanMinVarFreqSnp >varscan.snps.vcf  2>varscan_snps.err

    java -jar \$JARPATH mpileup2indel $resultPileup --vcf-sample-list vcf_sample_name --output-vcf 1 --p-value $params.varscanPValue --min-coverage $params.minCoverage --min-var-freq $params.varscanMinVarFreqCons >varscan.indels.vcf  2> varscan_indels.err

    java -jar \$JARPATH mpileup2cns $resultPileup --p-value $params.varscanPValue --min-coverage $params.minCoverage --min-var-freq $params.varscanMinVarFreqCons > varscan.cons 2>varscan_cons.err

    bgzip varscan.snps.vcf
    tabix -fp vcf varscan.snps.vcf.gz
    bgzip varscan.indels.vcf
    tabix -fp vcf varscan.indels.vcf.gz

    perl /usr/bin/parseVarscanToCoverage.pl \\
        --file varscan.cons \\
        --percentCutoff 60 \\
        --coverageCutoff $params.minCoverage \\
        --outputFile ${sampleName}.coverage.txt

    perl /usr/bin/maskGenome.pl \\
      -p $resultPileup \\
      -f $genomeReorderedFastaIndex \\
      -dc $params.minCoverage \\
      -o masked.fa

    fold -w 60 masked.fa > genome_masked.fa
    """

  stub:
    """
    touch varscan.snps.vcf.gz
    touch varscan.snps.vcf.gz.tbi
    touch varscan.indels.vcf.gz
    touch varscan.indels.vcf.gz.tbi
    touch genome_masked.fa
    touch varscan.cons.gz
    """
}

process concatSnpsAndIndels {
  container = 'biocontainers/bcftools:v1.9-1-deb_cv1'

  input:
    tuple val(sampleName), path(varscanSnpsVcfGz), path(varscanSnpsVcfGzTbi), path(varscanIndelsVcfGz), path(varscanIndelsVcfGzTbi), path(genomeMaskedFasta)

  output:
    tuple val(sampleName), path('varscan.concat.vcf'), path('genome_masked.fa')

  script:
    """
    set -euo pipefail
    bcftools concat \\
      -a \\
      -o varscan.concat.vcf $varscanSnpsVcfGz $varscanIndelsVcfGz
    """

  stub:
    """
    touch varscan.concat.vcf
    touch genome_masked.fa
    """

}

process makeCombinedVarscanIndex {
  container = 'veupathdb/dnaseqanalysis:1.0.0'

   publishDir "$params.outputDir", pattern: "*.concat.vcf.gz", mode: "copy"
   publishDir "$params.outputDir", pattern: "*.concat.vcf.gz.tbi", mode: "copy"

  input:
    tuple val(sampleName), path(varscanConcatVcf), path(genomeMaskedFasta)

  output:
    tuple val(sampleName), path('*.concat.vcf.gz'), path('*.concat.vcf.gz.tbi'), path('genome_masked.fa')

  script:
    """
    set -euo pipefail
    mv $varscanConcatVcf ${sampleName}.concat.vcf
    bgzip ${sampleName}.concat.vcf
    tabix -fp vcf ${sampleName}.concat.vcf.gz
    """

  stub:
    """
    touch varscan.concat.vcf.gz
    touch varscan.concat.vcf.gz.tbi
    touch genome_masked.fa
    """

}

process filterIndels {
  container = 'biocontainers/vcftools:v0.1.16-1-deb_cv1'

  input:
    tuple val(sampleName), path(varscanConcatVcfGz), path(varscanConcatVcfGzTbi), path(genomeMaskedFasta)

  output:
    tuple val(sampleName), path('output.recode.vcf')

  script:
    """
    set -euo pipefail
    vcftools \\
        --gzvcf $varscanConcatVcfGz \\
        --keep-only-indels \\
        --out output \\
        --recode
    """

  stub:
    """
    touch output.recode.vcf
    """

}

process makeIndelTSV {
  container = 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir", pattern: "output.tsv", mode: "copy", saveAs: { filename -> "${sampleName}.indel.tsv" }

  input:
    tuple val(sampleName), path(outputRecodeVcf)

  output:
    path('output.tsv')

  script:
    """
    set -euo pipefail
    findValues.pl \\
       -i $outputRecodeVcf \\
       -s ${sampleName} \\
       -o output.tsv
    """

  stub:
    """
    touch output.tsv
    """

}

process mergeVcfs {
  container = 'biocontainers/bcftools:v1.9-1-deb_cv1'

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

process makeMergedVarscanIndex {
  container = 'veupathdb/dnaseqanalysis:1.0.0'

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

process bcftoolsConsensus {
  container = 'biocontainers/bcftools:v1.9-1-deb_cv1'

  input:
    tuple val(sampleName), path(varscanConcatVcfGz), path(varscanConcatVcfGzTbi), path(genomeMaskedFasta)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    tuple val(sampleName), path('cons.fa')

  script:
    """
    set -euo pipefail
    bcftools consensus \\
      -I \\
      -f $genomeMaskedFasta $varscanConcatVcfGz > cons.fa
    """

  stub:
    """
    touch cons.fa
    """

}

process addSampleToDefline {
  container = 'veupathdb/dnaseqanalysis:1.0.0'

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
