#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process fastqc {
  container 'biocontainers/fastqc:v0.11.9_cv7'

  input:
    tuple val(sampleName), path(sampleFile)

  output:
    tuple val(sampleName), path('fastqc_output', type:'dir')

  script:
    """
    set -euo pipefail
    mkdir fastqc_output
    fastqc -o fastqc_output --extract $sampleFile
    """

  stub:
    """
    mkdir fastqc_output
    """
}

process fastqc_check {
  container 'veupathdb/shortreadaligner:1.0.0'

  input:
    tuple val(sampleName), path(sampleFile), path(fastqc_output)

  output:
    tuple val(sampleName), path('mateAEncoding')

  script:
    """
    set -euo pipefail
    fastqc_check.pl $fastqc_output mateAEncoding
    """

  stub:
    """
    touch mateAEncoding
    """

}

process trimmomatic {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple val(sampleName), path(sampleFile), path('mateAEncoding'), val(isPaired)

  output:
    tuple val(sampleName), path('sample_1P'), path('sample_2P')

  script:
    """
    set -euo pipefail

    if [ "$isPaired" = true ]; then
        mateAEncoding=\$(<mateAEncoding)

        if [ "$params.trimmomaticAdaptorsFile" = "NA" ]; then
            java -jar /usr/share/java/trimmomatic.jar PE \\
                -trimlog trimLog.txt $sampleFile -\$mateAEncoding \\
                -baseout sample ILLUMINACLIP:/usr/local/bin/All_adaptors-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
        else
            java -jar /usr/share/java/trimmomatic.jar PE \\
                -trimlog trimLog.txt $sampleFile -\$mateAEncoding \\
                -baseout sample ILLUMINACLIP:$params.trimmomaticAdaptorsFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
        fi
    else
        touch sample_2P
        mateAEncoding=\$(<mateAEncoding)

        if [ "$params.trimmomaticAdaptorsFile" = "NA" ]; then
            java -jar /usr/share/java/trimmomatic.jar SE \\
                -trimlog trimLog.txt $sampleFile \\
                -\$mateAEncoding sample_1P ILLUMINACLIP:/usr/local/bin/All_adaptors-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
        else
            java -jar /usr/share/java/trimmomatic.jar SE \\
                -trimlog trimLog.txt $sampleFile \\
                -\$mateAEncoding sample_1P ILLUMINACLIP:$params.trimmomaticAdaptorsFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
        fi
    fi
    """

  stub:
    """
    touch sample_1P
    touch sample_2P
    """

}
