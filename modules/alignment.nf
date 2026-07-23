#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bwaIndex {
  container 'veupathdb/dnaseqanalysis:1.0.1'

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
  container 'veupathdb/dnaseqanalysis:1.0.1'

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
          if zcat -f $sample_1p | normalizeReadNames.sh needs-normalizing \\
             || zcat -f $sample_2p | normalizeReadNames.sh needs-normalizing; then
              # Read names carry a mate suffix bwa-mem2 cannot auto-strip
              # (.1/.2 from SRA, :a/:b from legacy Illumina; bwa only handles /1,/2).
              # Strip them on the fly so mem_sam_pe does not abort with
              # "paired reads have different names". Stripping forces bwa to pair
              # positionally, so first guard that the mate files are in lockstep;
              # fail loud rather than silently mispair.
              c1=\$(zcat -f $sample_1p | wc -l)
              c2=\$(zcat -f $sample_2p | wc -l)
              if [ "\$c1" -ne "\$c2" ]; then
                  echo "ERROR: mate files for ${sampleName} differ in length (\$c1 vs \$c2 lines); cannot positionally pair reads." >&2
                  exit 1
              fi
              bwa-mem2 mem \\
                  -t $params.bwaThreads \\
                  -R '@RG\\tID:${sampleName}\\tSM:${sampleName}\\tPL:ILLUMINA' \\
                  genomeIndex \\
                  <(zcat -f $sample_1p | normalizeReadNames.sh) \\
                  <(zcat -f $sample_2p | normalizeReadNames.sh) \\
                  | samtools collate -o output.bam - && \
                  samtools fixmate -m output.bam fix.bam && \
                  samtools sort -o sort.bam fix.bam && \
                  samtools markdup -r sort.bam result_sorted.bam
          else
              # Read names already conform (identical, or standard /1,/2 that bwa
              # auto-strips); hand the raw files straight to bwa untouched.
              bwa-mem2 mem \\
                  -t $params.bwaThreads \\
                  -R '@RG\\tID:${sampleName}\\tSM:${sampleName}\\tPL:ILLUMINA' \\
                  genomeIndex \\
                  $sample_1p \\
                  $sample_2p \\
                  | samtools collate -o output.bam - && \
                  samtools fixmate -m output.bam fix.bam && \
                  samtools sort -o sort.bam fix.bam && \
                  samtools markdup -r sort.bam result_sorted.bam
          fi
      else
          bwa-mem2 mem \\
              -t $params.bwaThreads \\
              -R '@RG\\tID:${sampleName}\\tSM:${sampleName}\\tPL:ILLUMINA' \\
              genomeIndex \\
              $sample_1p \\
              | samtools sort -o result_sorted.bam -
      fi
      """

  stub:
    """
    touch result_sorted.bam
    """
}

process reorderFasta {
  container 'veupathdb/shortreadaligner:1.0.1'

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
    def jvmMem = (task.ext.memGb - 1) as int
    """
    set -euo pipefail
    JARPATH="/usr/picard/picard.jar"
    # GATK requires read group tags to be present in the BAM header
    java -Xmx${jvmMem}g -jar \$JARPATH AddOrReplaceReadGroups I=$resultSortedDsBam O=picard.bam RGID=$sampleName RGSM=$sampleName RGLB=NA RGPL=ILLUMINA RGPU=NA
    # GATK requires a sequence dictionary alongside the reference FASTA
    java -Xmx${jvmMem}g -jar \$JARPATH CreateSequenceDictionary R=$genomeReorderedFasta UR=$genomeReorderedFasta
    java -Xmx${jvmMem}g -jar \$JARPATH BuildBamIndex I=picard.bam
    java -Xmx${jvmMem}g -jar \$JARPATH CollectAlignmentSummaryMetrics R=$genomeReorderedFasta I=picard.bam O=summaryMetrics.txt
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
    def jvmMem = (task.ext.memGb - 1) as int
    """
    set -euo pipefail
    JARPATH="/usr/GenomeAnalysisTK.jar"
    # Scan the BAM to find genomic intervals where reads show signs of indel misalignment
    java -Xmx${jvmMem}g -jar \$JARPATH \\
      -I $picardBam \\
      -R $genomeReorderedFasta \\
      -T RealignerTargetCreator \\
      -allowPotentiallyMisencodedQuals \\
      -o forIndelRealigner.intervals 2>realaligner.err
    # Locally realign reads within the target intervals to produce a cleaner BAM
    java -Xmx${jvmMem}g -jar \$JARPATH \\
      -I $picardBam \\
      -R $genomeReorderedFasta \\
      -T IndelRealigner -targetIntervals forIndelRealigner.intervals \\
      -allowPotentiallyMisencodedQuals \\
      -o ${sampleName}.bam 2>indelRealigner.err

    mv ${sampleName}.bai ${sampleName}.bam.bai
    """

  stub:
    """
    touch ${sampleName}.bam
    touch ${sampleName}.bam.bai
    """
}

// Counts raw input fragments before trimming/filtering.
// For paired-end data, counts R1 only (one count per fragment).
process rawReadCount {
  container 'veupathdb/dnaseqanalysis:1.0.1'

  input:
    tuple val(sampleName), path(reads)

  output:
    tuple val(sampleName), path("${sampleName}.rawcount.tsv")

  script:
    def read1 = reads instanceof List ? reads[0] : reads
    """
    set -euo pipefail
    count=\$(( \$(zcat ${read1} | wc -l) / 4 ))
    printf 'sample\traw_fastq_reads\n' > ${sampleName}.rawcount.tsv
    printf '%s\t%s\n' "${sampleName}" "\$count" >> ${sampleName}.rawcount.tsv
    """

  stub:
    """
    touch ${sampleName}.rawcount.tsv
    """
}

// Extracts per-sample alignment metrics from the samtools stats SN (summary numbers) section:
//   trimmed_total_sequences, reads mapped, mapped read percentage, average read length
process samtoolsStats {
  container 'veupathdb/dnaseqanalysis:1.0.1'
    
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

    printf 'sample\ttrimmed_total_sequences\treads_mapped\tmapped_read_pct\taverage_read_length\n' \
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

// Joins the three per-sample stat files into a single published TSV for that sample.
// Columns: sample, raw_fastq_reads, trimmed_total_sequences, reads_mapped,
//          trimmed_mapped_pct, raw_mapped_pct, average_read_length, mean_coverage
process makeAlignmentStats {
  container 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  publishDir "${params.outputDir}/${sampleName}", mode: "copy"

  input:
    tuple val(sampleName), path(samtoolsFile), path(genomecovFile), path(rawcountFile)

  output:
    path "${sampleName}_alignment_stats.tsv"

  script:
    """
    set -euo pipefail

    # Each input file has a header line followed by exactly one data line.
    read -r _ trimmed_total reads_mapped trimmed_pct avg_length \
      < <(tail -n 1 $samtoolsFile)
    read -r _ mean_cov \
      < <(tail -n 1 $genomecovFile)
    read -r _ raw_reads \
      < <(tail -n 1 $rawcountFile)

    raw_pct=\$(awk -v rm="\$reads_mapped" -v rr="\$raw_reads" \
      'BEGIN { printf "%.4f", (rr > 0) ? rm / rr : 0 }')

    printf 'sample\traw_fastq_reads\ttrimmed_total_sequences\treads_mapped\ttrimmed_mapped_pct\traw_mapped_pct\taverage_read_length\tmean_coverage\n' \
      > ${sampleName}_alignment_stats.tsv
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${sampleName}" "\$raw_reads" "\$trimmed_total" "\$reads_mapped" \
      "\$trimmed_pct" "\$raw_pct" "\$avg_length" "\$mean_cov" \
      >> ${sampleName}_alignment_stats.tsv
    """

  stub:
    """
    touch ${sampleName}_alignment_stats.tsv
    """
}
