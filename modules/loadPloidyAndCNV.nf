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


workflow loadPloidyAndCNV {

  take:
    cnv_qch
    ploidy_qch
    
  main:
    loadPloidy(ploidy_qch)
    loadGeneCNV(cnv_qch)

}
