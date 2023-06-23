#!/usr/bin/env nextflow
nextflow.enable.dsl=2


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


process mergeVcfs {
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
    path 'variationFeature.dat', emit: variationFile
    path 'allele.dat', emit: alleleFile
    path 'product.dat', emit: productFile
  
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

process addFeatureIdsToVariation {
  container = 'veupathdb/dnaseqanalysis'
  
  publishDir "$params.outputDir", mode: "copy", pattern: 'variationFeatureFinal.dat'
  
  input:
    path variationFile
    path gusConfig
  
  output:
    path 'variationFeatureFinal.dat'
  
  script:
    template 'addFeatureIdsToVariation.bash'

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
    template 'insertVariation.bash'

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
    template 'insertProduct.bash'

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
    template 'insertAllele.bash'

  stub:
    """
    echo "insert allele"
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


workflow me {
 
  take:

    fastas_qch
    vcfs_qch
    indels_qch
    coverage_qch
    bw_qch
    bam_qch

  main:

    bigwigs = bw_qch.collectFile(storeDir: params.webServicesDir)
    bams = bam_qch.collectFile(storeDir: params.webServicesDir)

    coverages = coverage_qch.collectFile(storeDir: params.varscan_directory)

    combinedFastagz = fastas_qch.collectFile(name: 'CombinedFasta.fa.gz')
    combinedIndels = indels_qch.collectFile(name: 'indel.tsv') 

    checkUniqueIds(combinedFastagz) 

    allvcfs = vcfs_qch.collect()

    mergedVcf = mergeVcfs(allvcfs)    
  
    makeSnpFileResults = makeSnpFile(mergedVcf)
    
    processSeqVarsResults = processSeqVars(makeSnpFileResults.snpFile, params.cacheFile, params.undoneStrains, params.organism_abbrev, params.reference_strain, params.varscan_directory, params.genomeFastaFile, combinedFastagz, combinedIndels, params.gtfFile, coverages, bigwigs, bams)

    addFeatureIdsToVariationResults = addFeatureIdsToVariation(processSeqVarsResults.variationFile, params.gusConfig)
    
    insertVariation(params.extDbRlsSpec, addFeatureIdsToVariationResults)
    insertProduct(processSeqVarsResults.productFile)
    insertAllele(processSeqVarsResults.alleleFile)

    //snpEff(mergedVcfResults, params.gtfFile, params.genomeFastaFile)

}
