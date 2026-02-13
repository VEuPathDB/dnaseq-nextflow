#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process writeIndelConfigFile {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple val(sampleName), path(indelFile)

  output:
    tuple val(sampleName), path("${sampleName}_indelConfig.txt"), path(indelFile)

  script:
    """
    set -euo pipefail

    writeStudyConfig.pl \\
      --name $params.experimentName \\
      --file $indelFile \\
      --outputFile ${sampleName}_indelConfig.txt \\
      --sourceIdType NASequence \\
      --protocol Indel \\
      --profileSetName "${sampleName} (Indel)"
    """

  stub:
    """
    touch ${sampleName}_indelConfig.txt
    """
}

process writePloidyConfigFile {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple val(sampleName), path(ploidyFile)

  output:
    tuple val(sampleName), path("${sampleName}_ploidyConfig.txt"), path(ploidyFile)

  script:
    """
    set -euo pipefail

    writeStudyConfig.pl \\
      --name $params.experimentName \\
      --file $ploidyFile \\
      --outputFile ${sampleName}_ploidyConfig.txt \\
      --sourceIdType NASequence \\
      --protocol Ploidy \\
      --profileSetName "${sampleName} (CNV)"
    """

  stub:
    """
    touch ${sampleName}_ploidyConfig.txt
    """
}

process writeCNVConfigFile {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple val(sampleName), path(geneCNVFile)

  output:
    tuple val(sampleName), path("${sampleName}_geneCNVConfig.txt"), path(geneCNVFile)

  script:
    """
    set -euo pipefail

    writeStudyConfig.pl \\
      --outputFile ${sampleName}_geneCNVConfig.txt \\
      --file $geneCNVFile \\
      --name $params.experimentName \\
      --protocol geneCNV \\
      --sourceIdType gene \\
      --profileSetName "${sampleName} - GeneCNV"
    """

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
    """
    set -euo pipefail

    ga ApiCommonData::Load::Plugin::InsertStudyResults \\
      --inputDir '.' \\
      --configFile '$configFile' \\
      --extDbSpec '$extDbSpec' \\
      --studyName '$sampleName' \\
      --commit

    echo "DONE"
    """

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
    """
    set -euo pipefail

    ga ApiCommonData::Load::Plugin::InsertStudyResults \\
      --inputDir '.' \\
      --configFile '$configFile' \\
      --extDbSpec '$extDbSpec' \\
      --studyName '$sampleName' \\
      --commit

    echo "DONE"
    """

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
    """
    set -euo pipefail

    ga ApiCommonData::Load::Plugin::InsertStudyResults \\
      --inputDir '.' \\
      --configFile '$configFile' \\
      --extDbSpec '$extDbSpec' \\
      --studyName '$sampleName' \\
      --commit

    echo "DONE"
    """

  stub:
    """
    touch stdout
    """
}
