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


process writePloidyConfigFile {
  container 'veupathdb/dnaseqanalysis'

  publishDir "$params.outputDir"

  input:
    tuple val(sampleName), path(ploidyFile)

  output:
    path "${sampleName}_ploidyConfig.txt", emit: ploidyConfig
    path ploidyFile, emit: ploidyFile
    val(sampleName), emit: sampleName

  script:
    template 'writePloidyConfigFile.bash'

  stub:
    """
    touch ${sampleName}_ploidyConfig.txt
    touch $ploidyFile
    """

}


process writeCNVConfigFile {
  container 'veupathdb/dnaseqanalysis'

  publishDir "$params.outputDir"

  input:
    tuple val(sampleName), path(geneCNVFile)

  output:
    path "${sampleName}_geneCNVConfig.txt", emit: cnvConfig
    path geneCNVFile, emit: cnvFile
    val(sampleName), emit: sampleName

  script:
    template 'writeCNVConfigFile.bash'

  stub:
    """
    touch ${sampleName}_geneCNVConfig.txt
    """

}


process loadPloidy {
  tag "plugin"
  
  input:
    path 'ploidyFile'
    path 'configFile'
    val sampleName
    val extDbSpec

  script:
    template 'insertStudyResults.bash'

  stub:
    """
    touch stdout
    """
}


process loadGeneCNV {
  tag "plugin"
  
  input:
    path 'geneCNVFile'
    path 'configFile'
    val sampleName
    val extDbSpec

  script:
    template 'insertStudyResults.bash'

  stub:
    """
    touch stdout
    """

}


workflow lc {

  take:
    tpm_qch
    
  main:
  
    calcPloidyCNVResults = calculatePloidyAndGeneCNV(tpm_qch, params.footprintFile, params.gusConfig, params.ploidy, params.taxonId)
    
    writePloidyConfigFileResults = writePloidyConfigFile(calcPloidyCNVResults.ploidy) 
    writeCNVConfigFileResults = writeCNVConfigFile(calcPloidyCNVResults.geneCNV)
    
    loadPloidy(writePloidyConfigFileResults.ploidyFile, writePloidyConfigFileResults.ploidyConfig, writePloidyConfigFileResults.sampleName, params.extDbSpec)
    loadGeneCNV(writeCNVConfigFileResults.cnvFile, writeCNVConfigFileResults.cnvConfig, writeCNVConfigFileResults.sampleName, params.extDbSpec)

}
