#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process hisat2Index {
  container = 'veupathdb/shortreadaligner:1.0.0'

  input:
   path genomeFasta
   val fromBam
   val createIndex

  output:
   path 'genomeIndex*.ht2', emit: ht2_files
   val 'genomeIndex' , emit: genome_index_name

  script:
    """
    set -euo pipefail

    if [ "$fromBam" = true ]; then
        touch genomeIndex.1.ht2
        samtools faidx $genomeFasta
    elif [ "$createIndex" = true ]; then
        hisat2-build $genomeFasta genomeIndex
        samtools faidx $genomeFasta
    else
        TMP=$params.hisat2Index
        FILES=\$TMP*
        for f in \$FILES; do cp "\$f" "genomeIndex\${f#\$TMP}" ; done
        samtools faidx $genomeFasta
    fi
    """

  stub:
    """
    touch genomeIndex1.ht2
    """
}

process hisat2 {
    container = 'veupathdb/shortreadaligner:1.0.0'

    input:
      tuple val(sampleName), path(sampleFile), path('mateAEncoding'), path(sample_1p), path(sample_2p)
      val hisat2_index
      path 'genomeIndex.*.ht2'
      val fromBam
      val isPaired

    output:
      tuple val(sampleName), path('result_sorted.bam')

    script:
      """
      set -euo pipefail

      if [ "$fromBam" = true ]; then
          samtools view -bS $sampleFile | samtools sort - > result_sorted.bam
      elif [ "$isPaired" = true ]; then
          mateAEncoding=\$(<mateAEncoding)
          hisat2 --no-spliced-alignment \\
              -k 1 \\
              -p $params.hisat2Threads \\
              -q --\$mateAEncoding \\
              -x $hisat2_index \\
              -1 $sample_1p \\
              -2 $sample_2p  \\
              | samtools collate -o output.bam -
              samtools fixmate -m output.bam fix.bam
              samtools sort -o sort.bam fix.bam
              samtools markdup -r sort.bam result_sorted.bam
      else
          mateAEncoding=\$(<mateAEncoding)
          hisat2 --no-spliced-alignment \\
              -k 1 \\
              -p $params.hisat2Threads \\
              -q --\$mateAEncoding \\
              -x $hisat2_index \\
              -U $sample_1p \\
              | samtools collate -o output.bam -
              samtools fixmate -m output.bam fix.bam
              samtools sort -o sort.bam fix.bam
              samtools markdup -r sort.bam result_sorted.bam
      fi
      """

  stub:
    """
    touch result_sorted.bam
    """
}

process reorderFasta {
  container = 'veupathdb/shortreadaligner:1.0.0'

  input:
    tuple val(sampleName), path(resultSortedBam)
    path genomeFasta

  output:
    tuple path('genome_reordered.fa'), path('genome_reordered.fa.fai')

  script:
    """
    set -euo pipefail
    for seq in \$(samtools view -H $resultSortedBam | grep '^@SQ' | cut -f 2); do echo \${seq#*SN:}; done > regions.txt
    samtools faidx $genomeFasta \$(cat regions.txt) > genome_reordered.fa
    samtools faidx genome_reordered.fa
    """

  stub:
    """
    touch genome_reordered.fa
    touch genome_reordered.fa.fai
    """

}

process subsample {
  container = 'veupathdb/shortreadaligner:1.0.0'

  input:
    tuple val(sampleName), path(resultSortedBam)

  output:
    tuple val(sampleName), path('result_sorted_ds.bam')

  script:
    """
    set -euo pipefail
    samtools index result_sorted.bam

    # number of mapped reads is col 3 from idxstats
    frac=\$( samtools idxstats $resultSortedBam | awk 'BEGIN {total=0} {total += \$3} END {frac=$params.maxNumberOfReads/total;  if (frac > 1) {print 1} else {print frac}}' )

    # this will subsample fraction of mapped reads
    if awk "BEGIN {exit !(\$frac >= 1)}"
     then
       ln -s $resultSortedBam result_sorted_ds.bam
    else
       samtools view -b -s \$frac $resultSortedBam > result_sorted_ds.bam
    fi

    lineCount=\$(wc -l result_sorted_ds.bam)

    if [\$lineCount = 0]; then
        exit 1
    fi
    """

  stub:
    """
    touch result_sorted_ds.bam
    """
}

process picard {
  container = 'broadinstitute/picard:2.25.0'

  input:
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)
    tuple val(sampleName), path(resultSortedDsBam)

  output:
    tuple val(sampleName), path('genome_reordered.dict'), path('picard.bam'), path('picard.bai'), emit: bam_and_dict
    tuple val(sampleName), path('summaryMetrics.txt'), emit: metrics

  script:
    """
    set -euo pipefail
    JARPATH="/usr/picard/picard.jar"
    java -jar \$JARPATH AddOrReplaceReadGroups I=$resultSortedDsBam O=picard.bam RGID=$sampleName RGSM=$sampleName RGLB=NA RGPL=NA RGPU=NA
    java -jar \$JARPATH CreateSequenceDictionary R=$genomeReorderedFasta UR=$genomeReorderedFasta
    java -jar \$JARPATH BuildBamIndex I=picard.bam
    java -jar \$JARPATH CollectAlignmentSummaryMetrics R=$genomeReorderedFasta I=picard.bam O=summaryMetrics.txt
    """

  stub:
    """
    touch genome_reordered.dict
    touch picard.bam
    touch picard.bai
    touch summaryMetrics.txt
    """

}

process gatk {
  container = 'broadinstitute/gatk3:3.8-1'

  publishDir "$params.outputDir", pattern: "*.bam", mode: "copy"
  publishDir "$params.outputDir", pattern: "*.bai", mode: "copy"

  input:
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)
    tuple val(sampleName), path(genomeReorderedDict), path(picardBam), path(picardBamIndex)

  output:
    tuple val(sampleName), path("${sampleName}.bam"), path("${sampleName}.bam.bai")

  script:
    """
    set -euo pipefail
    JARPATH="/usr/GenomeAnalysisTK.jar"
    java -jar \$JARPATH \\
      -I $picardBam \\
      -R $genomeReorderedFasta \\
      -T RealignerTargetCreator \\
      -o forIndelRealigner.intervals 2>realaligner.err
    java -jar \$JARPATH \\
      -I $picardBam \\
      -R $genomeReorderedFasta \\
      -T IndelRealigner -targetIntervals forIndelRealigner.intervals \\
      -o ${sampleName}.bam 2>indelRealigner.err

    mv ${sampleName}.bai ${sampleName}.bam.bai
    """

  stub:
    """
    touch result_sorted_gatk.bam
    touch result_sorted_gatk.bai
    """
}
