#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process checkUniqueIds {

  input:
    path fastaDir

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
  publishDir "$params.outputDir", mode: "copy", pattern: 'merged.vcf.gz'

  input:
    path '*.vcf.gz'

  output:
    path 'merged.vcf.gz'
    path 'merge.vcf', emit: mergedVcf

  script:
    template 'mergeVcfs.bash'

  stub:
    """
    touch merged.vcf.gz
    touch merge.vcf
    """

}


process makeSnpFile {
  publishDir "$params.outputDir", mode: "copy"

  input:
    path 'merged.vcf'

  output: 
    path 'snpFile.tsv', emit: snpFile

  script:
    """
    perl /usr/bin/makeSnpFile.pl --vcf merged.vcf --output snpFile.tsv
    """

  stub:
    """
    touch snpFile.tsv
    """

}


process processSeqVars {
  publishDir "$params.outputDir", mode: "copy"

  input:
    path 'snpFile.tsv'
    path 'gusConfig.txt'
    path 'cache.txt'
    path 'undoneStrains.txt'
    val  transcript_extdb_spec
    val  organism_abbrev
    val  reference_strain
    val  extdb_spec
    path 'varscan_directory'
    path 'genome.fa'
    path 'consensus.fa.gz'
  
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
  publishDir "$params.outputDir", mode: "copy"

  input:
    path 'merged.vcf'
    path 'genes.gtf.gz'
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


workflow mergeExperiments {

  take:

    fastas_qch
    vcfs_qch
    
  main:
    
    checkResults = checkUniqueIds(params.fastaDir) 

    combinedFastagz = fastas_qch.collectFile(name: 'CombinedFasta.fa.gz', storeDir: params.outputDir )

    allvcfs = vcfs_qch.collect()

    mergeVcfsResults = mergeVcfs(allvcfs)
    
    makeSnpFileResults = makeSnpFile(mergeVcfsResults.mergedVcf)
    
    processSeqVars(makeSnpFileResults.snpFile, params.gusConfig, params.cacheFile, params.undoneStrains, params.transcript_extdb_spec, params.organism_abbrev, params.reference_strain, params.extdb_spec, params.varscan_directory, params.genomeFasta, combinedFastagz)

    snpEff(mergeVcfsResults.mergedVcf, params.databaseFile, params.sequenceFile)

}
