#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process writeIndelConfigFile {
  container = 'veupathdb/dnaseqanalysis'

  input:
    tuple val(sampleName), path(indelFile)

  output:
    tuple val(sampleName), path("${sampleName}_indelConfig.txt"), path(indelFile)
    
  script:
    template 'writeIndelConfigFile.bash'

  stub:
    """
    touch ${sampleName}_indelConfig.txt
    """
}


process writePloidyConfigFile {
  container = 'veupathdb/dnaseqanalysis'

  input:
    tuple val(sampleName), path(ploidyFile)

  output:
    tuple val(sampleName), path "${sampleName}_ploidyConfig.txt", path(ploidyFile)
    
  script:
    template 'writePloidyConfigFile.bash'

  stub:
    """
    touch ${sampleName}_ploidyConfig.txt
    """
}


process writeCNVConfigFile {
  container = 'veupathdb/dnaseqanalysis'

  input:
    tuple val(sampleName), path(geneCNVFile)

  output:
    tuple val(sampleName), path "${sampleName}_geneCNVConfig.txt", path(geneCNVFile)

  script:
    template 'writeCNVConfigFile.bash'

  stub:
    """
    touch ${sampleName}_geneCNVConfig.txt
    """
}


process loadIndels {
  tag "plugin"

  input:
    tuple val(sampleName), path(configFile), path(indelFile)
    val extDbSpec

  script:
    template 'insertStudyResults.bash'

  stub:
    """
    touch stdout
    """
}


process loadPloidy {
  tag "plugin"

  input:
    tuple val(sampleName), path(configFile), path(ploidyFile)
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
    tuple val(sampleName), path(configFile), path(geneCNVFile)
    val extDbSpec

  script:
    template 'insertStudyResults.bash'

  stub:
    """
    touch stdout
    """
}


workflow ls {

  take:
    indels_qch
    ploidy_qch
    cnv_qch
 
  main:
    writeIndelConfigFileResults = writeIndelConfigFile(indels_qch)
    writePloidyConfigFileResults = writePloidyConfigFile(ploidy_qch)
    writeCNVConfigFileResults = writeCNVConfigfile(cnv_qch)
    loadIndels(writeIndelConfigFileResults, params.extDbRlsSpec)
    loadPloidy(writePloidyConfigFileResults, params.extDbRlsSpec)
    loadGeneCNV(writeGeneCNVConfigFileResults, params.extDbRlsSpec)
}