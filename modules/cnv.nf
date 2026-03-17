#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process genomecov {
  container 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple val(sampleName), path(resultSortedGatkBam), path(resultSortedGatkIndex)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    tuple val(sampleName), path('coverage.bed')

  script:
    """
    set -euo pipefail
    bedtools genomecov \\
      -bg \\
      -ibam $resultSortedGatkBam \\
      -g $genomeReorderedFastaIndex >coverage.bed
    """

  stub:
    """
    touch coverage.bed
    """

}


process bedGraphToBigWig {
    container 'quay.io/biocontainers/ucsc-bedgraphtobigwig:469--h9b8f530_0'

    publishDir "$params.outputDir", mode: "copy"

    input:
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)
    tuple val(sampleName), path(coverageBed)

    output:
    path "${sampleName}.bw"

    script:
    """
    LC_COLLATE=C sort -k1,1 -k2,2n $coverageBed > sorted_input.bedgraph

    bedGraphToBigWig sorted_input.bedgraph $genomeReorderedFastaIndex ${sampleName}.bw
    """

    stub:
    """
    touch coverage.bw
    """
}

process normaliseCoverageToBigWig {
    container 'quay.io/biocontainers/ucsc-bedgraphtobigwig:469--h9b8f530_0'

    publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}_normalisedCoverage.bw" }

    input:
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)
    tuple val(sampleName), path(coverageBed)

    output:
    path "${sampleName}.bw"

    script:
    """
    LC_COLLATE=C sort -k1,1 -k2,2n $coverageBed > sorted_input.bedgraph

    bedGraphToBigWig sorted_input.bedgraph $genomeReorderedFastaIndex ${sampleName}.bw
    """

    stub:
    """
    touch ${sampleName}.bw
    """
}


process sortForCounting {
  container 'veupathdb/shortreadaligner:1.0.0'

    input:
    tuple val(sampleName), path(resultSortedGatkBam), path(resultSortedGatkBamIndex)

  output:
    tuple val(sampleName), path('result_sortByName.bam')

  script:
    """
    set -euo pipefail
    samtools sort -n $resultSortedGatkBam > result_sortByName.bam
    """

  stub:
    """
    touch result_sortByName.bam
    """

}

process htseqCount {
  container 'biocontainers/htseq:v0.11.2-1-deb-py3_cv1'

  //publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}.counts" }

  input:
    tuple val(sampleName), path(resultSortByNameBam)
    path gtfFile

  output:
    tuple val(sampleName), path('counts.txt')

  script:
    """
    set -euo pipefail
    htseq-count \\
      -f bam \\
      -s no \\
      -t CDS \\
      -i gene_id \\
      -a 0 $resultSortByNameBam $gtfFile > counts.txt
    """

  stub:
    """
    touch counts.txt
    """
}

process calculateTPM {
  container 'veupathdb/shortreadaligner:1.0.0'

  input:
    tuple val(sampleName), path(counts)
    path geneFootprintFile

  output:
    tuple val(sampleName), path('*.tpm')

  script:
    """
    set -euo pipefail
    makeTpmFromHtseqCountsCNV.pl \\
      --geneFootprintFile $geneFootprintFile \\
      --countFile $counts \\
      --tpmFile out.tpm
    """

  stub:
    """
    touch tpm.txt
    """

}

process makeWindowFile {
  container 'veupathdb/shortreadaligner:1.0.0'

  input:
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)
    val winLen

  output:
    tuple path('windows.bed'), path('genome.txt')

  script:
    """
    set -euo pipefail
    makeWindowedBed.pl \\
      --samtoolsIndex $genomeReorderedFastaIndex \\
      --winLen $winLen
    """

  stub:
    """
    touch windows.bed
    touch genome.txt
    """

}

process bedtoolsWindowed {
  container 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple path(windows), path(genome)
    tuple val(sampleName), path(resultSortedGatkBam), path(resultSortedGatkBamIndex)

  output:
    tuple val(sampleName), path('windowedCoverage.bed')

  script:
    """
    set -euo pipefail
    bedtools coverage \\
      -counts \\
      -sorted \\
      -g $genome \\
      -a $windows \\
      -b $resultSortedGatkBam > windowedCoverage.bed
    """

  stub:
    """
    touch windowedCoverage.bed
    """

}

process normaliseCoverage {
  container 'veupathdb/shortreadaligner:1.0.0'

  //publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}.bed" }

  input:
    tuple val(sampleName), path(windowedCoverage), path(summaryMetrics)

  output:
    tuple val(sampleName), path('normalisedCoverage.bed')

  script:
    """
    set -euo pipefail
    # NOTE final processing requires querying the DB so can stay in ReFlow
    normaliseCoverageCNV.pl \\
      --bedFile $windowedCoverage \\
      --summaryMetrics $summaryMetrics
    """

  stub:
    """
    touch normalisedCoverage.bed
    """
}

process makeSnpDensity {
  container  'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple val(sampleName), path(freebayesVcfGz), path(freebayesVcfGzTbi), path(snpsVcfGz), path(snpsVcfGzTbi), path(indelsVcfGz), path(indelsVcfGzTbi), path(gvcfGz), path(gvcfGzTbi)
    tuple path(windows), path(genome)

  output:
    tuple val(sampleName), path('snpDensity.bed'), path('indelDensity.bed')

  script:
    """
    set -euo pipefail
    zcat $snpsVcfGz | bedtools coverage \\
      -a $windows \\
      -b stdin -sorted \\
      -g $genome \\
      -counts > snpDensity.bed
    zcat $indelsVcfGz | bedtools coverage \\
      -a $windows \\
      -b stdin -sorted \\
      -g $genome \\
      -counts > indelDensity.bed
    """

  stub:
    """
    touch snpDensity.bed
    touch indelDensity.bed
    """
}

process makeDensityBigwigs {

  container 'quay.io/biocontainers/ucsc-bedgraphtobigwig:469--h9b8f530_0'

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}_${filename}" }

  input:
    tuple val(sampleName), path(snpDensity), path(indelDensity)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    tuple path('snpDensity.bw'), path('indelDensity.bw')

  script:
    """
    set -euo pipefail
    LC_COLLATE=C sort -k1,1 -k2,2n $snpDensity > sorted.snpDensity.bed
    LC_COLLATE=C sort -k1,1 -k2,2n $indelDensity > sorted.indelDensity.bed
    bedGraphToBigWig sorted.snpDensity.bed $genomeReorderedFastaIndex snpDensity.bw
    bedGraphToBigWig sorted.indelDensity.bed $genomeReorderedFastaIndex indelDensity.bw
    """

  stub:
    """
    touch snpDensity.bw
    touch indelDensity.bw
    """
}

process getHeterozygousSNPs {
  container 'veupathdb/vcf_parser_cnv:1.0.0'

  input:
    tuple val(sampleName), path(freebayesVcfGz), path(freebayesVcfGzTbi), path(snpsVcfGz), path(snpsVcfGzTbi), path(indelsVcfGz), path(indelsVcfGzTbi), path(gvcfGz), path(gvcfGzTbi)

  output:
    tuple val(sampleName), path('heterozygousSNPs.vcf')

  script:
    """
    set -euo pipefail
    makeHeterozygosityPlot.py --vcfFile $snpsVcfGz
    """

  stub:
    """
    touch heterozygousSNPs.vcf
    """
}

process makeHeterozygousDensityBed {
  container 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple path(windows), path(genome)
    tuple val(sampleName), path(heterozygousSNPs)

  output:
    tuple val(sampleName), path('heterozygousDensity.bed')

  script:
    """
    set -euo pipefail
    bedtools coverage \\
      -a $windows \\
      -b $heterozygousSNPs \\
      -sorted \\
      -g $genome \\
      -counts > heterozygousDensity.bed
    """

  stub:
    """
    touch heterozygousDensity.bed
    """
}

process makeHeterozygousDensityBigwig {
  container 'veupathdb/shortreadaligner:1.0.0'

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}_LOH.bw" }

  input:
    tuple val(sampleName), path(heterozygousDensityBed)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    path('heterozygousDensity.bw')

  script:
    """
    set -euo pipefail
    LC_COLLATE=C sort -k1,1 -k2,2n $heterozygousDensityBed > sorted.heterozygousDensity.bed
    bedGraphToBigWig sorted.heterozygousDensity.bed $genomeReorderedFastaIndex heterozygousDensity.bw
    """

  stub:
    """
    touch heterozygousDensity.bw
    """
}

process calculatePloidyAndGeneCNV {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> filename.endsWith("_CNVestimations.tsv") ? null : filename }

  input:
    tuple val(sampleName), path(sampleFile)
    path footprints
    val ploidy
    val taxonId
    path geneSourceIdOrtholog
    path chrsForCalc

  output:
    tuple val(sampleName), path( "${sampleName}_Ploidy.txt" ), emit: ploidy
    tuple val(sampleName), path( "${sampleName}_geneCNVs.txt" ), emit: geneCNV
    path "${sampleName}_CNVestimations.tsv"

  script:
    """
    set -euo pipefail

    calculatePloidy.pl \\
        --outputDir . \\
        --fpkmFile $sampleFile \\
        --sampleName $sampleName \\
        --taxonId  $taxonId \\
        --geneFootprints $footprints \\
        --ploidy $ploidy \\
        --chrsForCalcsFile $chrsForCalc

    calculateGeneCNVs.pl \\
        --fpkmFile $sampleFile \\
        --ploidy $ploidy \\
        --outputDir . \\
        --sampleName $sampleName \\
        --taxonId $taxonId \\
        --geneFootPrints $footprints \\
        --geneSourceIdOrthologFile $geneSourceIdOrtholog \\
        --chrsForCalcsFile $chrsForCalc
    """

  stub:
    """
    touch ${sampleName}_Ploidy.txt
    touch ${sampleName}_geneCNVs.txt
    touch ${sampleName}_CNVestimations.tsv
    """
}
