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

  stub:
    """
    touch stdout
    """
}

workflow ls {

  take:
    indels_qch
 
  main:
    loadIndels(indels_qch, params.extDbRlsSpec, params.genomeExtDbRlsSpec)

}