#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process copyBWToWS {
  tag "plugin"

  publishDir params.webServicesDir

  input:
    path(bw)

  output:
    path(bw)

  script:
    """
    """
}


process copyBAMToWS {
  tag "plugin"

  publishDir params.webServicesDir

  input:
    path(bam)

  output:
    path(bam)

  script:
    """
    """
}


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

    // TODO:  do we need a process to make the webservice directory?
    copyBWToWS(bw_qch)
    copyBAMToWS(bam_qch)
}
