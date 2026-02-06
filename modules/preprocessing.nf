#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process downloadBAMFromEBI {
  container = 'veupathdb/dnaseqanalysis:1.0.0'
  input:
    val id

  output:
    tuple val(id), path("${id}.bam")

  script:
    """
    set -euo pipefail
    wget --ftp-user $params.ebiFtpUser --ftp-password $params.ebiFtpPassword ftp://ftp-private.ebi.ac.uk:/upload/EBI/DNASeq/PlasmoDB_test/$params.organismAbbrev/Broad_HTS_Isolates_QuerySRA/$id/results.bam
    mv results.bam ${id}.bam
    """
}

process downloadFiles {
  container = 'veupathdb/humann:1.0.0'
  input:
    tuple val(strain), val(idList)

  output:
    tuple val(strain), path("${strain}**.fastq"), emit: files
    env isPaired, emit: isPaired

  script:
    """
    set -euo pipefail

    perl /usr/local/bin/getFilesFromSra.pl --strain $strain --idList $idList

    if [ -f "${strain}_2.fastq" ]; then
        export isPaired="true"
    else
        export isPaired="false"
    fi
    """
}

process fastqc {
  container = 'biocontainers/fastqc:v0.11.9_cv7'

  input:
    tuple val(sampleName), path(sampleFile)
    val fromBam

  output:
    tuple val(sampleName), path('fastqc_output', type:'dir')

  script:
    """
    set -euo pipefail

    if [ "$fromBam" = true ]; then
        mkdir fastqc_output
    else
        mkdir fastqc_output
        fastqc -o fastqc_output --extract $sampleFile
    fi
    """

  stub:
    """
    mkdir fastqc_output
    """
}

process fastqc_check {
  container = 'veupathdb/shortreadaligner:1.0.0'

  input:
    tuple val(sampleName), path(sampleFile), path(fastqc_output)
    val fromBam

  output:
    tuple val(sampleName), path('mateAEncoding')

  script:
    """
    set -euo pipefail

    if [ "$fromBam" = true ]; then
        touch mateAEncoding
    else
        fastqc_check.pl $fastqc_output mateAEncoding
    fi
    """

  stub:
    """
    touch mateAEncoding
    """

}

process trimmomatic {
  container = 'veupathdb/dnaseqanalysis:1.0.0'

  input:
    tuple val(sampleName), path(sampleFile), path('mateAEncoding')
    val fromBam
    val isPaired

  output:
    tuple val(sampleName), path('sample_1P'), path('sample_2P')

  script:
    """
    set -euo pipefail

    if [ "$fromBam" = true ]; then
        touch sample_1P
        touch sample_2P
    elif [ "$isPaired" = true ]; then
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
