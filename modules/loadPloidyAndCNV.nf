#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process calculatePloidyAndGeneCNV {
  publishDir "$params.outputDir"
  
  input:
    tuple val(sampleName), path(sampleFile)
    path 'footprints'
    path 'gusConfig'
    val ploidy
    val taxonId
  output:
    path '$sampleName_Ploidy.txt'
    path '$sampleName_geneCNVs.txt'
    path '$sampleName_CNVestimations.tsv'
  script:
    template 'calculatePloidyAndGeneCNV.bash'
}


process loadPloidy {
  tag "plugin"
  
  input:
    path(ploidy)

  script:
    template 'loadPloidy.bash'
}


process loadGeneCNV {
  tag "plugin"
  
  input:
    path(geneCNV)

  script:
    template 'loadGeneCNV.bash'
}


workflow loadPloidyAndCNV {

  take:
    tpm_qch
    
  main:
    calculatePloidyAndGeneCNV(tpm_qch, params.footprintFile, params.gusConfig, params.ploidy, params.taxonId)
    //loadPloidy()
    //loadGeneCNV()

}
