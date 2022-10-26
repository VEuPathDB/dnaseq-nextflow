#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process loadPloidy {
  tag "plugin"
  
  input:
    stdin

  output:
    stdout

  script:
    template 'loadPloidy.bash'
}


process loadGeneCNV {
  tag "plugin"
  
  input:
    stdin

  output:
    stdout

  script:
    template 'loadGeneCNV.bash'
}


process copyBWToWS {
  tag "plugin"
  
  input:
    stdin

  output:
    stdout

  script:
    template 'copyBWToWS.bash'
}


process copyBAMToWS {
  tag "plugin"

  input:
    stdin

  output:
    stdout

  script:
    template 'copyBAMToWS.bash'
}


process loadIndels {
  tag "plugin"

  input:
    stdin

  output:
    stdout

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
    copyBWToWS(bw_qch)
    copyBAMToWS(bam_qch)
    loadIndels(indels_qch)

}
