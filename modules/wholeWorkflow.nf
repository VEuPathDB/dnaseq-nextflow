#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//-------------------- ProcessSingleExperiment --------------------------

process downloadBAMFromEBI {
  container = 'veupathdb/dnaseqanalysis'
  input:
    val id

  output:
    tuple val(id), path("${id}.bam")         

  script:
    template 'downloadBAMFromEBI.bash'
}

process downloadFiles {
  container = 'veupathdb/humann'
  input:
    tuple val(strain), val(idList)

  output:
    tuple val(strain), path("${strain}**.fastq"), emit: files
    env isPaired, emit: isPaired

  script:
    template 'downloadFiles.bash'
}

process hisat2Index {
  container = 'veupathdb/shortreadaligner'

  input:
   path genomeFasta  
   val fromBam
   val createIndex

  output:
   path 'genomeIndex*.ht2', emit: ht2_files
   val 'genomeIndex' , emit: genome_index_name

  script:
    template 'hisat2Index.bash'

  stub:
    """
    touch genomeIndex1.ht2
    """
}

process fastqc {
  container = 'biocontainers/fastqc:v0.11.9_cv7'

  input:
    tuple val(sampleName), path(sampleFile)
    val fromBam

  output:
    tuple val(sampleName), path('fastqc_output', type:'dir')

  script:
      template 'fastqc.bash'

  stub:
    """
    mkdir fastqc_output
    """
}

process fastqc_check {
  container = 'veupathdb/shortreadaligner'

  input:
    tuple val(sampleName), path(sampleFile), path(fastqc_output)
    val fromBam

  output:
    tuple val(sampleName), path('mateAEncoding') 

  script:
    template 'fastqcCheck.bash'

  stub:
    """
    touch mateAEncoding
    """
}

process trimmomatic {
  container = 'veupathdb/shortreadaligner'

  input:
    tuple val(sampleName), path(sampleFile), path('mateAEncoding') 
    val fromBam
    val isPaired

  output:
    tuple val(sampleName), path('sample_1P'), path('sample_2P') 

  script:
    template 'trimmomatic.bash'

  stub:
    """
    touch sample_1P
    touch sample_2P
    """    
}

process hisat2 {
    container = 'veupathdb/shortreadaligner'

    input:
      tuple val(sampleName), path(sampleFile), path('mateAEncoding'), path(sample_1p), path(sample_2p) 
      val hisat2_index 
      path 'genomeIndex.*.ht2'
      val fromBam
      val isPaired

    output:
      tuple val(sampleName), path('result_sorted.bam')

    script:
      template 'hisat2.bash'
      
  stub:
    """
    touch result_sorted.bam
    """
}

process reorderFasta {
  container = 'veupathdb/shortreadaligner'
   
  input:
    tuple val(sampleName), path(resultSortedBam)
    path genomeFasta

  output:
    tuple path('genome_reordered.fa'), path('genome_reordered.fa.fai')

  script:
    template 'reorderFasta.bash'
     
  stub:
    """
    touch genome_reordered.fa
    touch genome_reordered.fa.fai
    """
}

process subsample {
  container = 'veupathdb/shortreadaligner'

  input:
    tuple val(sampleName), path(resultSortedBam)

  output:
    tuple val(sampleName), path('result_sorted_ds.bam')

  script:
    template 'subsample.bash'

  stub:
    """
    touch result_sorted_ds.bam
    """
}

process picard {
  container = 'broadinstitute/picard:2.25.0'

  input:
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)
    tuple val(sampleName), path(resultSortedDsBam) 

  output:
    tuple val(sampleName), path('genome_reordered.dict'), path('picard.bam'), path('picard.bai'), emit: bam_and_dict
    tuple val(sampleName), path('summaryMetrics.txt'), emit: metrics

  script:
    template 'picard.bash'

  stub:
    """
    touch genome_reordered.dict
    touch picard.bam
    touch picard.bai
    touch summaryMetrics.txt
    """
}

process gatk {
  container = 'broadinstitute/gatk3:3.8-1'

  publishDir "$params.outputDir", pattern: "*.bam", mode: "copy" 
  publishDir "$params.outputDir", pattern: "*.bai", mode: "copy" 

  input:
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)
    tuple val(sampleName), path(genomeReorderedDict), path(picardBam), path(picardBamIndex)

  output:
    tuple val(sampleName), path('${sampleName}.bam'), path('${sampleName}.bai'), emit: bamTuple
    path '${sampleName}.bam', emit: bamFile

  script:
    template 'gatk.bash'

  stub:
    """
    touch ${sampleName}.bam
    touch ${sampleName}.bai
    """
}

process mpileup {
  container = 'veupathdb/shortreadaligner'

  publishDir "$params.outputDir", pattern: "result.pileup", mode: "copy", saveAs: { filename -> "${sampleName}.result.pileup" }

  input:
    tuple val(sampleName), path (resultSortedGatkBam), path(resultSortedGatkBamIndex)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    tuple val(sampleName), path('result.pileup') 
    
  script:
    template 'mpileup.bash'

  stub:
    """
    touch result.pileup
    """
}

process varscan {
  container = 'veupathdb/dnaseqanalysis'

  publishDir "$params.outputDir/varscanCons", pattern: "${sampleName}.coverage.txt", mode: "copy" 

  input:
    tuple val(sampleName), path (resultSortedGatkBam), path(resultSortedGatkBamIndex), path(resultPileup) 
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi'), path('varscan.indels.vcf.gz'), path('varscan.indels.vcf.gz.tbi'), path('genome_masked.fa'), emit: vcf_files
    path '${sampleName}.coverage.txt', emit: coverageFile

  script:
    template 'varscan.bash'

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
    template 'concatSnpsAndIndels.bash'

  stub:
    """
    touch varscan.concat.vcf
    touch genome_masked.fa
    """
}

process makeCombinedVarscanIndex {
  container = 'veupathdb/dnaseqanalysis'
  
   publishDir "$params.outputDir", pattern: "*.concat.vcf.gz", mode: "copy"
   publishDir "$params.outputDir", pattern: "*.concat.vcf.gz.tbi", mode: "copy"

  input:
    tuple val(sampleName), path(varscanConcatVcf), path(genomeMaskedFasta) 

  output:
    tuple val(sampleName), path('*.concat.vcf.gz'), path('*.concat.vcf.gz.tbi'), path('genome_masked.fa')

  script:
    template 'makeCombinedVarscanIndex.bash'

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
    template 'filterIndels.bash'

  stub:
    """
    touch output.recode.vcf
    """
}

process makeIndelTSV {
  container = 'veupathdb/dnaseqanalysis'

  publishDir "$params.outputDir", pattern: "output.tsv", mode: "copy", saveAs: { filename -> "${sampleName}.indel.tsv" }

  input:
    tuple val(sampleName), path(outputRecodeVcf)

  output:
    path('output.tsv')

  script:
    template 'makeIndelTSV.bash'

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
    template 'mergeVcfsProcessSingle.bash'

  stub:
    """
    touch result.vcf.gz
    """
}

process makeMergedVarscanIndex {
  container = 'veupathdb/dnaseqanalysis'

  publishDir "$params.outputDir", mode: "copy"

  input:
    path(resultVcfGz) 

  output:
    tuple path('result.vcf.gz'), path('result.vcf.gz.tbi')

  script:
    template 'makeMergedVarscanIndex.bash'

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
    template 'bcftoolsConsensus.bash'

  stub:
    """
    touch cons.fa
    """
}

process addSampleToDefline {
  container = 'veupathdb/dnaseqanalysis'

  publishDir "$params.outputDir", mode: "copy", saveAs: { filename -> "${sampleName}_consensus.fa.gz" }

  input:
  tuple val(sampleName), path(consFasta)

  output:
    path 'unique_ids.fa.gz'

  script:
    template 'addSampleToDefline.bash'

  stub:
    """
    touch unique_ids.fa.gz
    """
}

process genomecov {
  container = 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple val(sampleName), path(resultSortedGatkBam), path(resultSortedGatkIndex) 
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    tuple val(sampleName), path('coverage.bed')

  script:
    template 'genomecov.bash'


  stub:
    """
    touch coverage.bed
    """
}

process bedGraphToBigWig {
  container = 'veupathdb/shortreadaligner'

  publishDir "$params.outputDir", mode: "copy" 

  input:
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)
    tuple val(sampleName), path(coverageBed) 

  output:
    path '${sampleName}.bw'
    
  script:
    template 'bedGraphToBigWig.bash'

  stub:
    """
    touch coverage.bw
    """
}

process sortForCounting {
  container = 'veupathdb/shortreadaligner'

    input:
    tuple val(sampleName), path(resultSortedGatkBam), path(resultSortedGatkBamIndex)

  output:
    tuple val(sampleName), path('result_sortByName.bam') 

  script:
    template 'sortForCounting.bash'

  stub:
    """
    touch result_sortByName.bam
    """
}

process htseqCount {
  container = 'biocontainers/htseq:v0.11.2-1-deb-py3_cv1'

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}.counts" }

  input:
    tuple val(sampleName), path(resultSortByNameBam) 
    path gtfFile 

  output: 
    tuple val(sampleName), path('counts.txt') 

  script:
    template 'htseqCount.bash'

  stub:
    """
    touch counts.txt
    """
}

process calculateTPM {
  container = 'veupathdb/shortreadaligner'

  publishDir "$params.outputDir/CNVs", mode: "copy"

  input:
    tuple val(sampleName), path(counts) 
    path geneFootprintFile 

  output:
    tuple val(sampleName), path('*.tpm')

  script:
    template 'calculateTpm.bash'

  stub:
    """
    touch tpm.txt
    """
}

process makeWindowFile {
  container = 'veupathdb/shortreadaligner'

  input:
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)
    val winLen 

  output:
    tuple path('windows.bed'), path('genome.txt') 

  script:
    template 'makeWindowFile.bash'

  stub:
    """
    touch windows.bed
    touch genome.txt
    """
}

process bedtoolsWindowed {
  container = 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple path(windows), path(genome) 
    tuple val(sampleName), path(resultSortedGatkBam), path(resultSortedGatkBamIndex) 

  output:
    tuple val(sampleName), path('windowedCoverage.bed') 

  script:
    template 'bedtoolsWindowed.bash'

  stub:
    """
    touch windowedCoverage.bed
    """
}

process normaliseCoverage {
  container = 'veupathdb/shortreadaligner' 

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}.bed" }

  input:
    tuple val(sampleName), path(windowedCoverage), path(summaryMetrics) 
   
  output:
    path 'normalisedCoverage.bed'

  script:
    template 'normaliseCoverage.bash'

  stub:
    """
    touch normalisedCoverage.bed
    """
}

process makeSnpDensity {
  container= 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple val(sampleName), path(varscanSnpsVcfGz), path(varscanSnpsVcfGzTbi), path(varscanIndelsVcfGz), path(varscanIndelsVcfGzTbi), path(genomeMaskedFasta) 
    tuple path(windows), path(genome) 

  output:
    tuple val(sampleName), path('snpDensity.bed'), path('indelDensity.bed')

  script: 
    template 'makeSnpDensity.bash'

  stub:
    """
    touch snpDensity.bed
    touch indelDensity.bed
    """
}

process makeDensityBigwigs {
  container = 'veupathdb/shortreadaligner'

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}_${filename}" }

  input:
    tuple val(sampleName), path(snpDensity), path(indelDensity)
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    tuple path('snpDensity.bw'), path('indelDensity.bw')

  script:
    template 'makeDensityBigWigs.bash'

  stub:
    """
    touch snpDensity.bw
    touch indelDensity.bw
    """
}

process getHeterozygousSNPs {
  container = 'veupathdb/vcf_parser_cnv'

  input:
    tuple val(sampleName), path(varscanSnpsVcfGz), path(varscanSnpsVcfGzTbi), path(varscanIndelsVcfGz), path(varscanIndelsVcfGzTbi), path(genomeMaskedFasta)

  output:
    tuple val(sampleName), path('heterozygousSNPs.vcf')

  script:
    template 'getHeterozygousSNPs.bash'

  stub:
    """
    touch heterozygousSNPs.vcf
    """
}

process makeHeterozygousDensityBed {
  container = 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple path(windows), path(genome) 
    tuple val(sampleName), path(heterozygousSNPs)

  output:
    tuple val(sampleName), path('heterozygousDensity.bed')

  script:
    template 'makeHeterozygousDensityBed.bash'

  stub:
    """
    touch heterozygousDensity.bed
    """
}

process makeHeterozygousDensityBigwig {
  container = 'veupathdb/shortreadaligner'

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}_LOH.bw" }

  input:
    tuple val(sampleName), path(heterozygousDensityBed) 
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)

  output:
    path('heterozygousDensity.bw') 

  script:
    template 'makeHeterozygousDensityBigWig.bash'


  stub:
    """
    touch heterozygousDensity.bw
    """
}

process calculatePloidyAndGeneCNV {
  container = 'veupathdb/dnaseqanalysis'

  publishDir "$params.outputDir", mode: "copy"
  
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
    template 'calculatePloidyAndGeneCNV.bash'

  stub:
    """
    touch ${sampleName}_Ploidy.txt
    touch ${sampleName}_geneCNVs.txt
    touch ${sampleName}_CNVestimations.tsv
    """
}

// ----------------- LoadSingleExperiment ------------------------------

process loadIndels {
  tag "plugin"

  input:
    tuple val(sampleName), path(indel)
    val extDbRlsSpec
    val genomeExtDbRlsSpec

  output:
    path "${sampleName}_indel.tsv"

  script:
    template 'loadIndelsUserDatasets.bash'

  stub:
    """
    touch stdout
    """
}

// ---------------- MergeExperiments -----------------------------------

process checkUniqueIds {
  container = 'veupathdb/dnaseqanalysis'
  
  input:
    path 'consensus.fa.gz'

  output:
    stdout

  script:
    template 'checkUniqueIds.bash'

  stub:
    """
    touch stdout
    """
}

process mergeVcfsMergeExperiment {
  container = 'veupathdb/dnaseqanalysis'
  publishDir "$params.outputDir", mode: "copy", pattern: 'merged.vcf.gz'

  input:
    path '*.vcf.gz'

  output:
    path 'merged.vcf.gz'

  script:
    template 'mergeVcfsMergeExperiments.bash'

  stub:
    """
    touch merged.vcf.gz
    touch merged.vcf
    """
}

process makeSnpFile {
  container = 'veupathdb/dnaseqanalysis'

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
  container = 'veupathdb/dnaseqanalysis'
  publishDir "$params.cacheFileDir", mode: "copy", pattern: "$params.cacheFile"
  publishDir "$params.outputDir", mode: "copy", pattern: 'allele.dat'
  publishDir "$params.outputDir", mode: "copy", pattern: 'product.dat'

  input:
    path snpFile
    path cacheFile
    path undoneStrainsFile
    val  organism_abbrev
    val  reference_strain
    path varscanDir
    path genomeFasta
    path consensusFasta
    path indelFile
    path gtfFile
    path coverageComplete
    path bigwigsComplete
    path bamsComplete
  
  output:
    path cacheFile
    path 'snpFeature.dat', emit: variationFile
    path 'allele.dat'
    path 'product.dat'
  
  script:
    template 'processSeqVars.bash'

  stub:
    """
    touch cache.txt
    touch snpFeature.dat
    touch allele.dat
    touch product.dat
    """
}

process addExtDbRlsIdToVariation {
  publishDir "$params.outputDir", mode: "copy"

  input:
    path variationFile
    val extDbSpec
    path gusConfig
  
  output:
    path 'variationFeature.dat'

  
  script:
    template 'addExtDbRlsId.bash'

  stub:
    """
    touch variationFeature.dat
    """
}

process snpEff {
  container = 'veupathdb/dnaseqanalysis'
  publishDir "$params.outputDir", mode: "copy"

  input:
    path 'merged.vcf'
    path 'genes.gtf'
    path 'sequences.fa.gz'

  output:
    path 'merged.ann.vcf'

  script:
    template 'snpEff.bash'    

  stub:
    """
    touch merged.ann.vcf
    """
}

// ----------------------- WORKFLOW ---------------------------

workflow wf {

  take:

    samples_qch

  main:

    genome_fasta_file = file(params.genomeFastaFile)

    hisat2IndexResults = hisat2Index(genome_fasta_file, params.fromBAM, params.createIndex)

    if(!params.local && !params.fromBAM) {

        downloadFilesResults = downloadFiles(samples_qch)

        fastqcResults = fastqc(downloadFilesResults.files, params.fromBAM)

        fastqc_checkResults = fastqc_check(downloadFilesResults.files.join(fastqcResults), params.fromBAM)

        trimmomaticResults = trimmomatic(downloadFilesResults.files.join(fastqc_checkResults), params.fromBAM, downloadFilesResults.isPaired)

        hisat2Results = hisat2(downloadFilesResults.files.join(fastqc_checkResults).join(trimmomaticResults), hisat2IndexResults.genome_index_name, hisat2IndexResults.ht2_files, params.fromBAM, downloadFilesResults.isPaired)
    }

    else if(!params.local && params.fromBAM) {

        files = downloadBAMFromEBI(samples_qch)

        fastqcResults = fastqc(files, params.fromBAM)

        fastqc_checkResults = fastqc_check(files.join(fastqcResults), params.fromBAM)

        trimmomaticResults = trimmomatic(files.join(fastqc_checkResults), params.fromBAM, 'NA')

        hisat2Results = hisat2(files.join(fastqc_checkResults).join(trimmomaticResults), hisat2IndexResults.genome_index_name, hisat2IndexResults.ht2_files, params.fromBAM, 'NA')
    
    }

    else {

        fastqcResults = fastqc(samples_qch, params.fromBAM)

        fastqc_checkResults = fastqc_check(samples_qch.join(fastqcResults), params.fromBAM)

        trimmomaticResults = trimmomatic(samples_qch.join(fastqc_checkResults), params.fromBAM, params.isPaired)

        hisat2Results = hisat2(samples_qch.join(fastqc_checkResults).join(trimmomaticResults), hisat2IndexResults.genome_index_name, hisat2IndexResults.ht2_files, params.fromBAM, params.isPaired)

    }
    
    reorderFastaResults = reorderFasta(hisat2Results.first(), genome_fasta_file)

    subsampleResults = subsample(hisat2Results)

    picardResults = picard(reorderFastaResults, subsampleResults)

    gatkResults = gatk(reorderFastaResults, picardResults.bam_and_dict )

    mpileupResults = mpileup(gatkResults.bamTuple, reorderFastaResults)

    varscanResults = varscan(gatkResults.bamTuple.join(mpileupResults), reorderFastaResults)
  
    concatSnpsAndIndelsResults = concatSnpsAndIndels(varscanResults.vcf_files)

    makeCombinedVarscanIndexResults = makeCombinedVarscanIndex(concatSnpsAndIndelsResults)

    filterIndelsResults = filterIndels(makeCombinedVarscanIndexResults)

    makeIndelTSVResults = makeIndelTSV(filterIndelsResults)

    // NOTE:  Must ensure the order here is consistent for the vcf files and their indexes;  the lists of paths are each sorted
    mergeVcfsResults = mergeVcfs(makeCombinedVarscanIndexResults.count(), makeCombinedVarscanIndexResults.map{ tuple it[1], it[2], "key" }.groupTuple(by: 2, sort: { a, b -> a <=> b } ))

    makeMergedVarscanIndexResults = makeMergedVarscanIndex(mergeVcfsResults)

    bcftoolsConsensusResults = bcftoolsConsensus(makeCombinedVarscanIndexResults, reorderFastaResults)

    addSampleToDeflineResults = addSampleToDefline(bcftoolsConsensusResults)

    genomecovResults = genomecov(gatkResults.bamTuple, reorderFastaResults)

    bedGraphToBigWigResults = bedGraphToBigWig(reorderFastaResults, genomecovResults)

    sortForCountingResults = sortForCounting(gatkResults.bamTuple)

    htseqCountResults = htseqCount(sortForCountingResults, params.gtfFile)

    calculateTPMResults = calculateTPM(htseqCountResults, params.footprintFile)

    calculatePloidyAndGeneCNV(calculateTPMResults, params.footprintFile, params.ploidy, params.taxonId, params.geneSourceIdOrthologFile, params.chrsForCalcFile)

    makeWindowFileResults = makeWindowFile(reorderFastaResults, params.winLen)

    bedtoolsWindowedResults =  bedtoolsWindowed(makeWindowFileResults, gatkResults.bamTuple)

    normaliseCoverageResults = normaliseCoverage(bedtoolsWindowedResults.join(picardResults.metrics))

    makeSnpDensityResults = makeSnpDensity(varscanResults.vcf_files, makeWindowFileResults)

    makeDensityBigwigsResults = makeDensityBigwigs(makeSnpDensityResults, reorderFastaResults)

    if (params.ploidy != 1) {

      getHeterozygousSNPsResults = getHeterozygousSNPs(varscanResults.vcf_files)

      makeHeterozygousDensityBedResults = makeHeterozygousDensityBed(makeWindowFileResults, getHeterozygousSNPsResults)

      makeHeterozygousDensityBigwig(makeHeterozygousDensityBedResults, reorderFastaResults)
    }

    loadIndelsResults = loadIndels(makeIndelTSVResults, params.extDbRlsSpec, params.genomeExtDbRlsSpec)

    bigwigs = bedGraphToBigWigResults.collectFile()
    bams = gatkResults.bamFile.collectFile()

    coverages = varscanResults.coverageFile.collectFile()

    combinedFastagz = addSampleToDeflineResults.collectFile(name: 'CombinedFasta.fa.gz')
    combinedIndels = loadIndelsResults.collectFile(name: 'indel.tsv')

    checkUniqueIds(combinedFastagz) 

    mergedVcf = mergeVcfsMergeExperiment(mergeVcfsResults)    
  
    makeSnpFileResults = makeSnpFile(mergedVcf)
    
    processSeqVarsResults = processSeqVars(makeSnpFileResults.snpFile, params.cacheFile, params.undoneStrains, params.organismAbbrev, params.reference_strain, "params.outputDir/varscanCons/", params.genomeFastaFile, combinedFastagz, combinedIndels, params.gtfFile, coverages, bigwigs, bams)

    addExtDbRlsIdToVariation(processSeqVarsResults.variationFile, params.extDbRlsSpec, params.gusConfig)

}
