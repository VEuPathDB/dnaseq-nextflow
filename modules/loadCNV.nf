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
    path "${sampleName}_ploidyConfig.txt", emit: ploidyConfig
    path ploidyFile, emit: ploidyFile

  script:
    template 'writePloidyConfigFile.bash'
}


process writeCNVConfigFile {
  container 'veupathdb/dnaseqanalysis'

  publishDir "$params.outputDir"

  input:
    tuple val(sampleName), path(geneCNVFile)

  output:
    path "${sampleName}_geneCNVConfig.txt", emit: cnvConfig
    path geneCNVFile, emit: cnvFile

  script:
    template 'writeCNVConfigFile.bash'
}


process loadPloidy {
  tag "plugin"
  
  input:
    path 'ploidyFile'
    path 'configFile'

  script:
    template 'insertStudyResults.bash'
}


process loadGeneCNV {
  tag "plugin"
  
  input:
    path 'geneCNVFile'
    path 'configFile'

  script:
    template 'insertStudyResults.bash'
}


workflow loadCNV {

  take:
    tpm_qch
    
  main:
  
    calcPloidyCNVResults = calculatePloidyAndGeneCNV(tpm_qch, params.footprintFile, params.gusConfig, params.ploidy, params.taxonId)
    
    writePloidyConfigFileResults = writePloidyConfigFile(calcPloidyCNVResults.ploidy) 
    writeCNVConfigFileResults = writeCNVConfigFile(calcPloidyCNVResults.geneCNV)
    
    loadPloidy(writePloidyConfigFileResults.ploidyFile, writePloidyConfigFileResults.ploidyConfig)
    loadGeneCNV(writeCNVConfigFileResults.cnvFile, writeCNVConfigFileResults.cnvConfig)

}
