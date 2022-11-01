#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process calculatePloidyAndGeneCNV {
  container 'veupathdb/dnaseqanalysis'

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
}

process writePloidyConfigFile {
  container 'veupathdb/dnaseqanalysis'
  publishDir "$params.outputDir"
  input:
    tuple val(sampleName), path(ploidyFile)
  output:
    path "${sampleName}_ploidyConfig.txt"
  script:
    template 'writePloidyConfigFile.bash'
}

process writeCNVConfigFile {
  container 'veupathdb/dnaseqanalysis'
  publishDir "$params.outputDir"
  input:
    tuple val(sampleName), path(geneCNVFile)
  output:
    path "${sampleName}_geneCNVConfig.txt"
  script:
    template 'writeCNVConfigFile.bash'
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
    calcPloidyCNVResults = calculatePloidyAndGeneCNV(tpm_qch, params.footprintFile, params.gusConfig, params.ploidy, params.taxonId)
    writePloidyConfigFile(calcPloidyCNVResults.ploidy)
    writeCNVConfigFile(calcPloidyCNVResults.geneCNV)
    //loadPloidy()
    //loadGeneCNV()

}
