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
    path 'merged.vcf', emit: mergedVcf

  script:
    template 'mergeVcfs.bash'

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
  publishDir "$params.outputDir", mode: "copy"

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
  
  output:
    path 'cache.txt'
    path 'snpFeature.dat'
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

  main:
    coverages = coverage_qch.collectFile(storeDir: params.varscan_directory)

    combinedFastagz = fastas_qch.collectFile(name: 'CombinedFasta.fa.gz')
    combinedIndels = indels_qch.collectFile(name: 'indel.tsv')

    checkUniqueIds(combinedFastagz) 

    mergedVcf = vcfs_qch.collect()

    //mergeVcfsResults = mergeVcfs(allvcfs)
    
    makeSnpFileResults = makeSnpFile(mergedVcf)
    
    processSeqVars(makeSnpFileResults.snpFile, params.cacheFile, params.undoneStrains, params.organism_abbrev, params.reference_strain, params.varscan_directory, params.genomeFastaFile, combinedFastagz, combinedIndels, params.gtfFile, coverages)

    //snpEff(mergedVcf, params.gtfFile, params.genomeFastaFile)

}
