#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow ls {

  take:
    
    bam_qch
    bw_qch
    
  main:

    bw_qch | collectFile(storeDir: params.webServicesDir)
    bam_qch | collectFile(storeDir: params.webServicesDir)
    
}
