#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process loadIndels {
  tag "plugin"

  input:
    path(indel)
    val extDbRlsSpec
    val genomeExtDbRlsSpec

  script:
    template 'loadIndels.bash'
}


workflow loadSingleExperiment {

  take:
    indels_qch
    bam_qch
    bw_qch
    
  main:

    loadIndels(indels_qch, params.extDbRlsSpec, params.genomeExtDbRlsSpec)

    bw_qch | collectFile(storeDir: params.webServicesDir)
    bam_qch | collectFile(storeDir: params.webServicesDir)
    
}
