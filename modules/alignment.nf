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
    bwa-mem2 index genomeIndex
    samtools faidx $genomeFasta
    """

  stub:
    """
    touch genomeIndex.amb genomeIndex.ann genomeIndex.bwt.2bit.64 genomeIndex.pac genomeIndex.0123
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
          bwa-mem2 mem \\
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
          bwa-mem2 mem \\
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
    touch ${sampleName}.bam
    touch ${sampleName}.bam.bai
    """
}

// Extracts per-sample alignment metrics from the samtools stats SN (summary numbers) section:
//   raw total sequences, reads mapped, mapped read percentage, average read length
process samtoolsStats {
  container 'veupathdb/dnaseqanalysis:1.0.0'
    
  input:
    tuple val(sampleName), path(bamFile), path(bamIndex)

  output:
    tuple val(sampleName), path("${sampleName}.samtools.tsv")

  script:
    """
    set -euo pipefail

    # samtools stats SN lines have the form:  SN  <label>:  <value>  [# comment]
    # grep ^SN filters to summary numbers; cut -f 3 takes the numeric value column
    samtools stats ${bamFile} > samtools_stats.txt
    raw_total=\$(grep '^SN.*raw total sequences:' samtools_stats.txt | cut -f 3)
    reads_mapped=\$(grep '^SN.*reads mapped:'      samtools_stats.txt | cut -f 3)
    avg_length=\$(grep   '^SN.*average length:'    samtools_stats.txt | cut -f 3)

    mapped_pct=\$(echo "\$reads_mapped \$raw_total" | awk '{printf "%.4f", \$1/\$2}')

    printf 'sample\traw_total_sequences\treads_mapped\tmapped_read_pct\taverage_read_length\n' \
      > ${sampleName}.samtools.tsv
    printf '%s\t%s\t%s\t%s\t%s\n' \
      "${sampleName}" "\$raw_total" "\$reads_mapped" "\$mapped_pct" "\$avg_length" \
      >> ${sampleName}.samtools.tsv
    """

  stub:
    """
    touch ${sampleName}.samtools.tsv
    """
}

// Computes mean genome coverage per sample using a bedtools genomecov frequency histogram.
// Without -bg, bedtools emits rows: chrom, depth, numBases, chromSize, fraction.
// The "genome" rows aggregate across all chromosomes.
// mean_coverage = sum(depth * numBases) / sum(numBases)  i.e. total bases covered / genome size
process bedtoolsGenomecovStats {
  container 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple val(sampleName), path(bamFile), path(bamIndex)
    tuple path(genomeFasta), path(genomeFastaIndex)

  output:
    tuple val(sampleName), path("${sampleName}.genomecov.tsv")

  script:
    """
    set -euo pipefail

    mean_cov=\$(bedtools genomecov -ibam ${bamFile} -g ${genomeFastaIndex} \
      | awk '
          /^genome/ {
            total_base_coverage += \$2 * \$3   # depth * numBases_at_that_depth
            genome_size         += \$3
          }
          END { printf "%.2f", total_base_coverage / genome_size }
        ')

    printf 'sample\tmean_coverage\n' > ${sampleName}.genomecov.tsv
    printf '%s\t%s\n' "${sampleName}" "\$mean_cov" >> ${sampleName}.genomecov.tsv
    """

  stub:
    """
    touch ${sampleName}.genomecov.tsv
    """
}

// Joins per-sample samtools and genomecov stats and publishes a single TSV for all samples.
// Columns: sample, raw_total_sequences, reads_mapped, mapped_read_pct, average_read_length, mean_coverage
process mergeAlignmentStats {
  container 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  publishDir "$params.outputDir", mode: "copy"

  input:
    path samtoolsFiles
    path genomecovFiles

  output:
    path 'alignment_stats.tsv'

  script:
    """
    set -euo pipefail

    # Join the two per-sample TSVs by sample name.
    # FILENAME matching routes columns into named arrays; END block prints joined rows.
    awk -F'\t' '
      FNR == 1 { next }   # skip header in each input file
      FILENAME ~ /samtools/ {
        raw_total[   \$1] = \$2
        reads_mapped[\$1] = \$3
        mapped_pct[  \$1] = \$4
        avg_length[  \$1] = \$5
      }
      FILENAME ~ /genomecov/ { mean_cov[\$1] = \$2 }
      END {
        for (sample in raw_total)
          print sample "\t" raw_total[sample] "\t" reads_mapped[sample] "\t" \
                mapped_pct[sample] "\t" avg_length[sample] "\t" mean_cov[sample]
      }
    ' *.samtools.tsv *.genomecov.tsv | sort > data.tsv

    printf 'sample\traw_total_sequences\treads_mapped\tmapped_read_pct\taverage_read_length\tmean_coverage\n' \
      > alignment_stats.tsv
    cat data.tsv >> alignment_stats.tsv
    """

  stub:
    """
    touch alignment_stats.tsv
    """
}
