#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bwaIndex {
  container 'veupathdb/dnaseqanalysis:1.0.0'

  input:
   path genomeFasta

  output:
   path 'genomeIndex.*', emit: index_files
   path genomeFasta, emit: genome_fasta

  script:
    """
    set -euo pipefail
    cp $genomeFasta genomeIndex
    bwa index genomeIndex
    samtools faidx $genomeFasta
    """

  stub:
    """
    touch genomeIndex.amb genomeIndex.ann genomeIndex.bwt genomeIndex.pac genomeIndex.sa
    """
}

process bwaMem {
  container 'veupathdb/dnaseqanalysis:1.0.0'

    input:
      tuple val(sampleName), path(sampleFile), path('mateAEncoding'), path(sample_1p), path(sample_2p), val(isPaired)
      path indexFiles
      path genomeFasta

    output:
      tuple val(sampleName), path('result_sorted.bam')

    script:
      """
      set -euo pipefail

      if [ "$isPaired" = true ]; then
          bwa mem \\
              -t $params.bwaThreads \\
              -R '@RG\\tID:${sampleName}\\tSM:${sampleName}\\tPL:ILLUMINA' \\
              genomeIndex \\
              $sample_1p \\
              $sample_2p \\
              | samtools collate -o output.bam -
              samtools fixmate -m output.bam fix.bam
              samtools sort -o sort.bam fix.bam
              samtools markdup -r sort.bam result_sorted.bam
      else
          bwa mem \\
              -t $params.bwaThreads \\
              -R '@RG\\tID:${sampleName}\\tSM:${sampleName}\\tPL:ILLUMINA' \\
              genomeIndex \\
              $sample_1p \\
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
  container 'veupathdb/shortreadaligner:1.0.0'

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


// Prepares the BAM file for GATK processing:
//   1. AddOrReplaceReadGroups   - adds required read group metadata (sample name, platform, etc.)
//   2. CreateSequenceDictionary - builds a .dict file from the reference FASTA that GATK requires
//   3. BuildBamIndex            - indexes the BAM for random access by genomic position
//   4. CollectAlignmentSummaryMetrics - generates QC stats (% mapped reads, mismatch rate, etc.)
process picard {
  container 'broadinstitute/picard:2.25.0'

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
    # GATK requires read group tags to be present in the BAM header
    java -jar \$JARPATH AddOrReplaceReadGroups I=$resultSortedDsBam O=picard.bam RGID=$sampleName RGSM=$sampleName RGLB=NA RGPL=NA RGPU=NA
    # GATK requires a sequence dictionary alongside the reference FASTA
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

// Performs local indel realignment using GATK3 to reduce false-positive SNP calls.
// BWA-MEM aligns reads independently, so reads spanning an indel can be positioned
// inconsistently. Realignment considers all reads in a region together to fix this.
//   1. RealignerTargetCreator - identifies intervals likely to harbor misaligned reads near indels
//   2. IndelRealigner         - locally realigns reads within those intervals
process gatk {
  container 'broadinstitute/gatk3:3.8-1'

  publishDir "$params.outputDir", pattern: "*.bam", mode: "copy"
  publishDir "$params.outputDir", pattern: "*.bai", mode: "copy"

  input:
    tuple path(genomeReorderedFasta), path(genomeReorderedFastaIndex)
    tuple val(sampleName), path(genomeReorderedDict), path(picardBam), path(picardBamIndex)

  output:
    tuple val(sampleName), path("${sampleName}.bam"), path("${sampleName}.bam.bai"), emit: bamTuple
    path("${sampleName}.bam"), emit: bamFiles

  script:
    """
    set -euo pipefail
    JARPATH="/usr/GenomeAnalysisTK.jar"
    # Scan the BAM to find genomic intervals where reads show signs of indel misalignment
    java -jar \$JARPATH \\
      -I $picardBam \\
      -R $genomeReorderedFasta \\
      -T RealignerTargetCreator \\
      -o forIndelRealigner.intervals 2>realaligner.err
    # Locally realign reads within the target intervals to produce a cleaner BAM
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
