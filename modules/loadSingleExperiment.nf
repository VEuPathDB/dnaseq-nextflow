#!/usr/bin/env nextflow
nextflow.enable.dsl=2


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


process copyBWToWS {
  tag "plugin"

    // TODO add a publishDir  for webservices here

  input:
    path(bw)

  output:
    path(bw)

  // TODO: may not need the script here.  just the dir to publish ?
  script:
    template 'copyBWToWS.bash'
}


process copyBAMToWS {
    tag "plugin"

    // TODO add a publishDir  for webservices here

    input:
    path(bam)

    output:
    path(bam)

    // TODO: may not need the script here.  just the dir to publish ?
    script:
    template 'copyBAMToWS.bash'
}


process loadIndels {
  tag "plugin"

  input:
    path(indel)

  script:
    template 'loadIndels.bash'
}


workflow loadSingleExperiment {

  take:
    indels_qch
    bam_qch
    bw_qch
    cnv_qch
    ploidy_qch
    
  main:
    loadPloidy(ploidy_qch)
    loadGeneCNV(cnv_qch)

    loadIndels(indels_qch)

    // TODO:  do we need a process to make the webservice directory?
    copyBWToWS(bw_qch)
    copyBAMToWS(bam_qch)
}