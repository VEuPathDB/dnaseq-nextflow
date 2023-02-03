#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process calculatePloidyAndGeneCNV {

  publishDir "$params.outputDir", pattern: "*_CNVestimations.tsv"
  
  input:
    tuple val(sampleName), path(sampleFile)
    path 'footprints'
    path 'gusConfig'
    val ploidy
    val taxonId

  output:
    tuple val(sampleName), path( "${sampleName}_Ploidy.txt" ), emit: ploidy
    tuple val(sampleName), path( "${sampleName}_geneCNVs.txt" ), emit: geneCNV
    path "${sampleName}_CNVestimations.tsv"

  script:
    template 'calculatePloidyAndGeneCNV.bash'

  stub:
    """
    touch ${sampleName}_Ploidy.txt
    touch ${sampleName}_geneCNVs.txt
    touch ${sampleName}_CNVestimations.tsv
    """
}


workflow lc {

  take:
    tpm_qch
    
  main:
  
    calcPloidyCNVResults = calculatePloidyAndGeneCNV(tpm_qch, params.footprintFile, params.gusConfig, params.ploidy, params.taxonId)
    
}
