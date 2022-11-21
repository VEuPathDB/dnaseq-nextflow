#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process hisat2Index {
  container = 'veupathdb/shortreadaligner'

  input:
   path 'genome.fa'  
    
  output:
   path 'genomeIndex*.ht2', emit: ht2_files
   val 'genomeIndex' , emit: genome_index_name

  script: 
    if(params.fromBAM)
      template 'hisat2IndexFromBam.bash'
    else if( params.createIndex )
      template 'hisat2IndexCreateIndex.bash'
    else
      template 'hisat2IndexNoBamNoIndex.bash'
}


process fastqc {
  container = 'biocontainers/fastqc:v0.11.9_cv7'

  input:
    tuple val(sampleName), path(sampleFile)

  output:
    tuple val(sampleName), path('fastqc_output', type:'dir')

  script:
    if(params.fromBAM)
      template 'fastqcBam.bash'
    else 
      template 'fastqcNoBam.bash'
}


process fastqc_check {
  container = 'veupathdb/shortreadaligner'

  input:
    tuple val(sampleName), path(sampleFile), path('fastqc_output') 

  output:
    tuple val(sampleName), path('mateAEncoding') 

  script:
    if(params.fromBAM)
      template 'fastqc_checkBam.bash'
    else 
      template 'fastqc_checkNoBam.bash'            
}


process trimmomatic {
  container = 'veupathdb/shortreadaligner'

  input:
    tuple val(sampleName), path(sampleFile), path('mateAEncoding') 
    path adaptorsFile 

  output:
    tuple val(sampleName), path('sample_1P'), path('sample_2P') 

  script:
    if(params.fromBAM)
      template 'trimmomaticBam.bash'      
    else if( params.isPaired )
      template 'trimmomaticPaired.bash'
    else 
      template 'trimmomaticSingle.bash'
}


process hisat2 {
    container = 'veupathdb/shortreadaligner'

    input:
      tuple val(sampleName), path(sampleFile), path('mateAEncoding'), path('sample_1p'), path('sample_2p') 
      val hisat2_index 
      path 'genomeIndex.*.ht2' 

    output:
      tuple val(sampleName), path('result_sorted.bam')

    script:
      if(params.fromBAM)
        template 'hisat2Bam.bash'
      else if( params.isPaired )
        template 'hisat2Paired.bash'
      else 
        template 'hisat2Single.bash'
}


process reorderFasta {
  container = 'veupathdb/shortreadaligner'
   
  input:
    tuple val(sampleName), path('result_sorted.bam')
    path 'genome.fa'

  output:
    tuple path('genome_reordered.fa'), path('genome_reordered.fa.fai')

  script:
    template 'reorderFasta.bash'
}


process subsample {
  container = 'veupathdb/shortreadaligner'

  input:
    tuple val(sampleName), path('result_sorted.bam')

  output:
    tuple val(sampleName), path('result_sorted_ds.bam')

  script:
    template 'subsample.bash'
}


process picard {
  container = 'broadinstitute/picard:2.25.0'

  input:
    tuple path('genome_reordered.fa'), path('genome_reordered.fa.fai')
    path picardJar
    tuple val(sampleName), path('result_sorted_ds.bam') 

  output:
    tuple val(sampleName), path('genome_reordered.dict'), path('picard.bam'), path('picard.bai'), emit: bam_and_dict
    tuple val(sampleName), path('summaryMetrics.txt'), emit: metrics

  script:
    template 'picard.bash'
}


process gatk {
  container = 'broadinstitute/gatk3:3.8-1'

  publishDir "$params.outputDir", pattern: "*.bam", mode: "copy", saveAs: { filename -> "${sampleName}.bam" }
  publishDir "$params.outputDir", pattern: "*.bai", mode: "copy", saveAs: { filename -> "${sampleName}.bam.bai" }

  input:
    path gatkJar 
    tuple path('genome_reordered.fa'), path('genome_reordered.fa.fai')
    tuple val(sampleName), path('genome_reordered.dict'), path('picard.bam'), path('picard.bai')

  output:
    tuple val(sampleName), path('result_sorted_gatk.bam'), path('result_sorted_gatk.bai')

  script:
    template 'gatk.bash'
}


process mpileup {
  container = 'veupathdb/shortreadaligner'

  publishDir "$params.outputDir", pattern: "result.pileup", mode: "copy", saveAs: { filename -> "${sampleName}.result.pileup" }

  input:
    tuple val(sampleName), path ('result_sorted_gatk.bam'), path('result_sorted_gatk.bam.bai')
    tuple path('genome_reordered.fa'), path('genome_reordered.fa.fai')

  output:
    tuple val(sampleName), path('result.pileup') 
    
  script:
    template 'mpileup.bash'
}


process varscan {
  container = 'veupathdb/dnaseqanalysis'

  publishDir "$params.outputDir", pattern: "varscan.cons.gz", mode: "copy", saveAs: { filename -> "${sampleName}.varscan.cons.gz" }

  input:
    tuple val(sampleName), path ('result_sorted_gatk.bam'), path('result_sorted_gatk.bam.bai'), path('result.pileup') 
    path varscanJar
    tuple path('genome_reordered.fa'), path('genome_reordered.fa.fai')

  output:
    tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi'), path('varscan.indels.vcf.gz'), path('varscan.indels.vcf.gz.tbi'), path('genome_masked.fa'), emit: vcf_files
    path 'varscan.cons.gz'

  script:
    template 'varscan.bash'
}


process concatSnpsAndIndels {
  container = 'biocontainers/bcftools:v1.9-1-deb_cv1'

  input:
    tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi'), path('varscan.indels.vcf.gz'), path('varscan.indels.vcf.gz.tbi'), path('genome_masked.fa')

  output:
    tuple val(sampleName), path('varscan.concat.vcf'), path('genome_masked.fa')

  script:
    template 'concatSnpsAndIndels.bash'
}


process makeCombinedVarscanIndex {
  container = 'veupathdb/dnaseqanalysis'
  
   publishDir "$params.outputDir", pattern: "varscan.concat.vcf.gz", mode: "copy", saveAs: { filename -> "${sampleName}.concat.vcf.gz" }
   publishDir "$params.outputDir", pattern: "varscan.concat.vcf.gz.tbi", mode: "copy", saveAs: { filename -> "${sampleName}.concat.vcf.gz.tbi" }

  input:
    tuple val(sampleName), path('varscan.concat.vcf'), path('genome_masked.fa') 

  output:
    tuple val(sampleName), path('varscan.concat.vcf.gz'), path('varscan.concat.vcf.gz.tbi'), path('genome_masked.fa')

  script:
    template 'makeCombinedVarscanIndex.bash'
}


process filterIndels {
  container = 'biocontainers/vcftools:v0.1.16-1-deb_cv1'

  input:
    tuple val(sampleName), path('varscan.concat.vcf.gz'), path('varscan.concat.vcf.gz.tbi'), path('genome_masked.fa')

  output:
    tuple val(sampleName), path('output.recode.vcf')

  script:
    template 'filterIndels.bash'
}


process makeIndelTSV {
  container = 'veupathdb/dnaseqanalysis'

  publishDir "$params.outputDir", pattern: "output.tsv", mode: "copy", saveAs: { filename -> "${sampleName}.indel.tsv" }

  input:
    tuple val(sampleName), path('output.recode.vcf')

  output:
    path('output.tsv')

  script:
    template 'makeIndelTSV.bash'
}


process mergeVcfs {
  container = 'biocontainers/bcftools:v1.9-1-deb_cv1'

  input:
    val vcfCount
    tuple path('*.vcf.gz'), path('*.vcf.gz.tbi'), val(key)

  output:
    path('result.vcf.gz') 
    
  script:
    if (vcfCount > 1)
      template 'mergeVcfsCount.bash'
    else
      template 'mergeVcfsNoCount.bash'
}


process makeMergedVarscanIndex {
  container = 'veupathdb/dnaseqanalysis'

  publishDir "$params.outputDir", mode: "copy"

  input:
    path('result.vcf.gz') 

  output:
    tuple path('result.vcf.gz'), path('result.vcf.gz.tbi')

  script:
    template 'makeMergedVarscanIndex.bash'
}


process bcftoolsConsensus {
  container = 'biocontainers/bcftools:v1.9-1-deb_cv1'

  input:
    tuple val(sampleName), path('varscan.concat.vcf.gz'), path('varscan.concat.vcf.gz.tbi'), path('genome_masked.fa')
    tuple path('genome_reordered.fa'), path('genome_reordered.fa.fai')

  output:
    tuple val(sampleName), path('cons.fa')
    
  script:
    template 'bcftoolsConsensus.bash'
}


process addSampleToDefline {
  container = 'veupathdb/dnaseqanalysis'

  publishDir "$params.outputDir", mode: "copy", saveAs: { filename -> "${sampleName}_consensus.fa.gz" }

  input:
  tuple val(sampleName), path('cons.fa')

  output:
    path 'unique_ids.fa.gz'

  script:
    template 'addSampleToDefline.bash'
}


process genomecov {
  container = 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple val(sampleName), path('result_sorted_gatk.bam'), path('result_sorted_gatk.bai') 
    tuple path('genome_reordered.fa'), path('genome_reordered.fa.fai')

  output:
    tuple val(sampleName), path('coverage.bed')

  script:
    template 'genomecov.bash'
}


process bedGraphToBigWig {
  container = 'veupathdb/shortreadaligner'

  publishDir "$params.outputDir", mode: "copy", saveAs: { filename -> "${sampleName}.bw" }

  input:
    tuple path('genome_reordered.fa'), path('genome_reordered.fa.fai')
    tuple val(sampleName), path('coverage.bed') 

  output:
    path 'coverage.bw'
    
  script:
    template 'bedGraphToBigWig.bash'
}


process sortForCounting {
  container = 'veupathdb/shortreadaligner'

    input:
    tuple val(sampleName), path('result_sorted_gatk.bam'), path('result_sorted_gatk.bam.bai')

  output:
    tuple val(sampleName), path('result_sortByName.bam') 

  script:
    template 'sortForCounting.bash'
}


process htseqCount {
  container = 'biocontainers/htseq:v0.11.2-1-deb-py3_cv1'

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}.counts" }

  input:
    tuple val(sampleName), path('result_sortByName.bam') 
    path 'gtfFile' 

  output: 
    tuple val(sampleName), path('counts.txt') 

  script:
    template 'htseqCount.bash'
}


process calculateTPM {
  container = 'veupathdb/shortreadaligner'

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}.tpm" }

  input:
    tuple val(sampleName), path('counts.txt') 
    path 'geneFootprintFile' 

  output:
    path('tpm.txt')

  script:
    template 'calculateTpm.bash'
}


process makeWindowFile {
  container = 'veupathdb/shortreadaligner'

  input:
    tuple path('genome_reordered.fa'), path('genome_reordered.fa.fai')
    val winLen 

  output:
    tuple path('windows.bed'), path('genome.txt') 

  script:
    template 'makeWindowFile.bash'
}


process bedtoolsWindowed {
  container = 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple path('windows.bed'), path('genome.txt') 
    tuple val(sampleName), path('result_sorted_gatk.bam'), path('result_sorted_gatk.bam.bai') 

  output:
    tuple val(sampleName), path('windowedCoverage.bed') 

  script:
    template 'bedtoolsWindowed.bash'
}


process normaliseCoverage {
  container = 'veupathdb/shortreadaligner' 

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}.bed" }

  input:
    tuple val(sampleName), path('windowedCoverage.bed'), path('summaryMetrics.txt') 
   
  output:
    path 'normalisedCoverage.bed'

  script:
    template 'normaliseCoverage.bash'
} 


process makeSnpDensity {
  container= 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi'), path('varscan.indels.vcf.gz'), path('varscan.indels.vcf.gz.tbi'), path('genome_masked.fa') 
    tuple path('windows.bed'), path('genome.txt') 

  output:
    tuple val(sampleName), path('snpDensity.bed'), path('indelDensity.bed')

  script: 
    template 'makeSnpDensity.bash'
}


process makeDensityBigwigs {
  container = 'veupathdb/shortreadaligner'

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}_${filename}" }

  input:
    tuple val(sampleName), path('snpDensity.bed'), path('indelDensity.bed')
    tuple path('genome_reordered.fa'), path('genome_reordered.fa.fai')

  output:
    tuple path('snpDensity.bw'), path('indelDensity.bw')

  script:
    template 'makeDensityBigWigs.bash'
}


process getHeterozygousSNPs {
  container = 'veupathdb/vcf_parser_cnv'

  input:
    tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi') 

  output:
    tuple val(sampleName), path('heterozygousSNPs.vcf') 

  script:
    template 'getHeterozygousSNPs.bash'
}


process makeHeterozygousDensityBed {
  container = 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple path('windows.bed'), path('genome.txt') 
    tuple val(sampleName), path('heterozygousSNPs.vcf')

  output:
    tuple val(sampleName), path('heterozygousDensity.bed')

  script:
    template 'makeHeterozygousDensityBed.bash'
}


process makeHeterozygousDensityBigwig {
  container = 'veupathdb/shortreadaligner'

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}_LOH.bw" }

  input:
    tuple val(sampleName), path('heterozygousDensity.bed') 
    tuple path('genome_reordered.fa'), path('genome_reordered.fa.fai')

  output:
    path('heterozygousDensity.bw') 

  script:
    template 'makeHeterozygousDensityBigWig.bash'
}


workflow processSingleExperiment {

  take:

    samples_qch

  main:
 
    gatk_jar_path = file(params.gatkJar)
    picard_jar_path = file(params.picardJar)
    varscan_jar_path = file(params.varscanJar)

    genome_fasta_file = file(params.genomeFastaFile)

    hisat2IndexResults = hisat2Index(genome_fasta_file)

    fastqcResults = fastqc(samples_qch)

    fastqc_checkResults = fastqc_check(samples_qch.join(fastqcResults))

    trimmomaticResults = trimmomatic(samples_qch.join(fastqc_checkResults), params.trimmomaticAdaptorsFile)

    hisat2Results = hisat2(samples_qch.join(fastqc_checkResults).join(trimmomaticResults), hisat2IndexResults.genome_index_name, hisat2IndexResults.ht2_files)

    reorderFastaResults = reorderFasta(hisat2Results.first(), genome_fasta_file)

    subsampleResults = subsample(hisat2Results)

    picardResults = picard(reorderFastaResults, picard_jar_path, subsampleResults)

    gatkResults = gatk(gatk_jar_path, reorderFastaResults, picardResults.bam_and_dict )

    mpileupResults = mpileup(gatkResults, reorderFastaResults)

    varscanResults = varscan(gatkResults.join(mpileupResults), varscan_jar_path, reorderFastaResults)
  
    concatSnpsAndIndelsResults = concatSnpsAndIndels(varscanResults.vcf_files)

    makeCombinedVarscanIndexResults = makeCombinedVarscanIndex(concatSnpsAndIndelsResults)

    filterIndelsResults = filterIndels(makeCombinedVarscanIndexResults)

    makeIndelTSV(filterIndelsResults)

    // NOTE:  Must ensure the order here is consistent for the vcf files and their indexes;  the lists of paths are each sorted
    mergeVcfsResults = mergeVcfs(makeCombinedVarscanIndexResults.count(), makeCombinedVarscanIndexResults.map{ tuple it[1], it[2], "key" }.groupTuple(by: 2, sort: { a, b -> a <=> b } ))

    makeMergedVarscanIndexResults = makeMergedVarscanIndex(mergeVcfsResults)

    bcftoolsConsensusResults = bcftoolsConsensus(makeCombinedVarscanIndexResults, reorderFastaResults)

    addSampleToDefline(bcftoolsConsensusResults)

    genomecovResults = genomecov(gatkResults, reorderFastaResults)

    bedgraphToBigWigResults = bedGraphToBigWig(reorderFastaResults, genomecovResults)

    sortForCountingResults = sortForCounting(gatkResults)

    htseqCountResults = htseqCount(sortForCountingResults, params.gtfFile)

    calculateTPMResults = calculateTPM(htseqCountResults, params.geneFootprintFile)

    makeWindowFileResults = makeWindowFile(reorderFastaResults, params.winLen)

    bedtoolsWindowedResults =  bedtoolsWindowed(makeWindowFileResults, gatkResults)

    normaliseCoverageResults = normaliseCoverage(bedtoolsWindowedResults.join(picardResults.metrics))

    makeSnpDensityResults = makeSnpDensity(varscanResults.vcf_files, makeWindowFileResults)

    makeDensityBigwigsResults = makeDensityBigwigs(makeSnpDensityResults, reorderFastaResults)

    if (params.ploidy != 1) {

      getHeterozygousSNPsResults = getHeterozygousSNPs(varscanResults.vcf_files)

      makeHeterozygousDensityBedResults = makeHeterozygousDensityBed(makeWindowFileResults, getHeterozygousSNPsResults)

      makeHeterozygousDensityBigwig(makeHeterozygousDensityBedResults, reorderFastaResults)
    }

}
