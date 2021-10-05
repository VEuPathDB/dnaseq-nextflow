if(params.fromBAM) {
    samples_qch = Channel.fromPath([params.inputDir + '/**/*.bam']).map { file -> tuple(file.baseName, [file]) }
}
else if(params.isPaired) {
   samples_qch = Channel.fromFilePairs([params.inputDir + '/**/*_{1,2}.fastq', params.inputDir + '/**/*_{1,2}.fastq.gz', params.inputDir + '/**/*_{1,2}.fq.gz'])
}
else {
    samples_qch = Channel.fromPath([params.inputDir + '/**/*.fastq', params.inputDir + '/**/*.fastq.gz', params.inputDir + '/**/*.fq.gz']).map { file -> tuple(file.baseName, [file]) }
}

samples_qch.into { samples_to_fastqc_qch; samples_to_fastqc_check_qch; samples_to_trimmomatic_qch; samples_to_hisat2_qch }

gatk_jar_path_ch = file(params.gatkJar);
picard_jar_path_ch = file(params.picardJar);
varscan_jar_path_ch = file(params.varscanJar);


process hisat2Index {
     container = 'veupathdb/shortreadaligner'

     input:
     	path 'genome.fa' from params.genomeFastaFile
    
    output:
	    path 'genomeIndex*.ht2' into hisat2_index_path_ch
        val 'genomeIndex' into hisat2_index_value_ch
        path 'genome.fa.fai' into genome_index_path_ch
        path 'genome.fa' into genome_fasta_path_ch

    script: 
    if(params.fromBAM)
      """
      touch genomeIndex.1.ht2
      samtools faidx genome.fa
      """
    else if( params.createIndex )
      """
      hisat2-build genome.fa genomeIndex
      samtools faidx genome.fa
      """
    else
      """
      TMP=$params.hisat2Index
      FILES=\$TMP*
      for f in \$FILES; do cp "\$f" "genomeIndex\${f#\$TMP}" ; done
      samtools faidx genome.fa
      """
}


process fastqc {
  container = 'biocontainers/fastqc:v0.11.9_cv7'

  input:
	tuple val(sampleName), path(sampleFile) from samples_to_fastqc_qch

  output:
        tuple val(sampleName), path('fastqc_output', type:'dir') into fastqc_dir_qch

   script:
    if(params.fromBAM)
	    """
        mkdir fastqc_output   
        """

    else 
	    """
        mkdir fastqc_output   
        fastqc -o fastqc_output --extract $sampleFile
        """
}


process fastqc_check {
  container = 'veupathdb/shortreadaligner'

  input:
	tuple val(sampleName), path(sampleFile), path('fastqc_output') from samples_to_fastqc_check_qch.join(fastqc_dir_qch)

  output:
        tuple val(sampleName), path('mateAEncoding') into (encoding_to_trimmomatic_qch, encoding_to_hisat2_qch)

  script:
    if(params.fromBAM)
       """
       touch mateAEncoding
       """

    else 
       """
       fastqc_check.pl fastqc_output mateAEncoding
       """
}


process trimmomatic {
  container = 'veupathdb/shortreadaligner'

  input:
	tuple val(sampleName), path(sampleFile), path('mateAEncoding') from samples_to_trimmomatic_qch.join(encoding_to_trimmomatic_qch)
    path adaptorsFile from params.trimmomaticAdaptorsFile

  output:
    tuple val(sampleName), path('sample_1P'), path('sample_2P') into trimmomatic_qch

  script:
    if(params.fromBAM)
        """
        touch sample_1P
        touch sample_2P
        """
    else if( params.isPaired )
        """
         mateAEncoding=\$(<mateAEncoding)
         java org.usadellab.trimmomatic.TrimmomaticPE -trimlog trimLog.txt $sampleFile -\$mateAEncoding -baseout sample ILLUMINACLIP:$adaptorsFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
        """
    else 
	    """
        touch sample_2P
        mateAEncoding=\$(<mateAEncoding)
        java org.usadellab.trimmomatic.TrimmomaticSE -trimlog trimLog.txt $sampleFile -\$mateAEncoding sample_1P ILLUMINACLIP:$adaptorsFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
        """
}


process hisat2 {
    container = 'veupathdb/shortreadaligner'

    input:
	    tuple val(sampleName), path(sampleFile), path('mateAEncoding'), path('sample_1p'), path('sample_2p') from samples_to_hisat2_qch.join(encoding_to_hisat2_qch).join(trimmomatic_qch)
        val hisat2_index from hisat2_index_value_ch
        path 'genomeIndex.*.ht2' from hisat2_index_path_ch

    output:
        tuple val(sampleName), path('result_sorted.bam') into hisat2_to_subsample_qch
    
    script:
    if(params.fromBAM)
        """
        samtools view -bS $sampleFile | samtools sort - > result_sorted.bam
        """

    else if( params.isPaired )
        """
        mateAEncoding=\$(<mateAEncoding)
        hisat2 --no-spliced-alignment -k 1 -p $params.hisat2Threads -q --\$mateAEncoding -x $hisat2_index -1 sample_1p -2 sample_2p  | samtools view -bS | samtools sort | samtools rmdup - result_sorted.bam

        """

    else 
	    """
        mateAEncoding=\$(<mateAEncoding)
        hisat2 --no-spliced-alignment -k 1 -p $params.hisat2Threads -q --\$mateAEncoding -x $hisat2_index -U sample_1p | samtools view -bS - | samtools sort - > result_sorted.bam
        """
}



process subsample {
    container = 'veupathdb/shortreadaligner'
    input:
       tuple val(sampleName), path('result_sorted.bam') from hisat2_to_subsample_qch

    output:
        tuple val(sampleName), path('result_sorted_ds.bam') into subsample_to_picard_qch

    script:
        """
        samtools index result_sorted.bam

        # number of mapped reads is col 3 from idxstats
        frac=\$( samtools idxstats result_sorted.bam | awk 'BEGIN {total=0} {total += \$3} END {frac=$params.maxNumberOfReads/total;  if (frac > 1) {print 1} else {print frac}}' )

        # this will subsample fraction of mapped reads
        if awk "BEGIN {exit !(\$frac >= 1)}"
        then
          ln -s result_sorted.bam result_sorted_ds.bam

        else
          samtools view -b -s \$frac result_sorted.bam > result_sorted_ds.bam
        fi
        """
}

process picard {
     container = 'broadinstitute/picard:2.25.0'

     input:
        path 'genome.fa' from genome_fasta_path_ch
        path 'genome.fa.fai' from genome_index_path_ch
        path picardJar from picard_jar_path_ch
        tuple val(sampleName), path('result_sorted_ds.bam') from subsample_to_picard_qch


     output:
        tuple val(sampleName), path('genome.dict'), path('picard_reordered.bam'), path('picard_reordered.bai') into picard_to_gatk_qch
        tuple val(sampleName), path('summaryMetrics.txt') into picard_to_normaliseCoverage_qch

    script:
        def jarPath = picardJar.name == 'NA' ? "/usr/picard/picard.jar" : picardJar.name
        """
        java -jar $jarPath AddOrReplaceReadGroups I=result_sorted_ds.bam O=picard.bam RGID=$sampleName RGSM=$sampleName RGLB=NA RGPL=NA RGPU=NA

        java -jar $jarPath CreateSequenceDictionary R=genome.fa UR=genome.fa
        java -jar $jarPath ReorderSam I=picard.bam O=picard_reordered.bam SD=genome.dict
        java -jar $jarPath BuildBamIndex I=picard_reordered.bam
        java -jar $jarPath CollectAlignmentSummaryMetrics R=genome.fa I=picard_reordered.bam O=summaryMetrics.txt
        """
}


process gatk {
    container = 'broadinstitute/gatk3:3.8-1'

    publishDir "$params.outputDir", pattern: "*.bam", mode: "copy", saveAs: { filename -> "${sampleName}.bam" }
    publishDir "$params.outputDir", pattern: "*.bai", mode: "copy", saveAs: { filename -> "${sampleName}.bam.bai" }

    input:
        path gatkJar from gatk_jar_path_ch
        path 'genome.fa' from genome_fasta_path_ch
        path 'genome.fa.fai' from genome_index_path_ch
        tuple val(sampleName), path('genome.dict'), path('picard_reordered.bam'), path('picard_reordered.bai') from picard_to_gatk_qch

    output:
        tuple val(sampleName), path('result_sorted_gatk.bam'), path('result_sorted_gatk.bai') into (gatk_bam_to_mpileup_qch, gatk_bam_to_varscan_qch, gatk_bam_to_genomecov_qch, gatk_to_sortForCounting_qch, gatk_to_bedtoolsWindowed_qch)

    script:
        def jarPath = gatkJar.name == 'NA' ? "/usr/GenomeAnalysisTK.jar" : gatkJar.name
        """
        java -jar $jarPath -I picard_reordered.bam -R genome.fa -T RealignerTargetCreator -o forIndelRealigner.intervals 2>realaligner.err
        java -jar $jarPath -I picard_reordered.bam -R genome.fa -T IndelRealigner -targetIntervals forIndelRealigner.intervals -o result_sorted_gatk.bam 2>indelRealigner.err
        """
}


process mpileup {
    container = 'veupathdb/shortreadaligner'

    input:
        tuple val(sampleName), path ('result_sorted_gatk.bam'), path('result_sorted_gatk.bam.bai') from gatk_bam_to_mpileup_qch
        path 'genome.fa' from genome_fasta_path_ch
        path 'genome.fa.fai' from genome_index_path_ch

    
    output:
        tuple val(sampleName), path('result.pileup') into mpileup_qch
    
    script:
	    """
        samtools mpileup -A -f genome.fa -B result_sorted_gatk.bam > result.pileup 2>pileup.err
        """
}


process varscan {
    container = 'veupathdb/dnaseqanalysis'

    publishDir "$params.outputDir", pattern: "varscan.cons.gz", mode: "copy", saveAs: { filename -> "${sampleName}.varscan.cons.gz" }

    input:
        tuple val(sampleName), path ('result_sorted_gatk.bam'), path('result_sorted_gatk.bam.bai'), path('result.pileup') from gatk_bam_to_varscan_qch.join(mpileup_qch)
        path varscanJar from varscan_jar_path_ch

    output:
        tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi'), path('varscan.indels.vcf.gz'), path('varscan.indels.vcf.gz.tbi') into (varscan_to_concat_qch, varscan_to_snpDensity_qch)
        tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi') into varscan_to_LOH_qch
        path('varscan.cons.gz')
    
    script:
        def jarPath = varscanJar.name == 'NA' ? "/usr/local/VarScan.jar" : varscanJar.name
	    """
       echo $sampleName >vcf_sample_name
       java -jar $jarPath mpileup2snp result.pileup --vcf-sample-list vcf_sample_name --output-vcf 1 --p-value $params.varscanPValue --min-coverage $params.varscanMinCoverage --min-var-freq $params.varscanMinVarFreqSnp >varscan.snps.vcf  2>varscan_snps.err
       java -jar $jarPath mpileup2indel result.pileup --vcf-sample-list vcf_sample_name --output-vcf 1 --p-value $params.varscanPValue --min-coverage $params.varscanMinCoverage --min-var-freq $params.varscanMinVarFreqCons >varscan.indels.vcf  2> varscan_indels.err
       java -jar $jarPath mpileup2cns result.pileup --p-value $params.varscanPValue --min-coverage $params.varscanMinCoverage --min-var-freq $params.varscanMinVarFreqCons > varscan.cons 2>varscan_cons.err
       bgzip varscan.snps.vcf
       tabix -fp vcf varscan.snps.vcf.gz
       bgzip varscan.indels.vcf
       tabix -fp vcf varscan.indels.vcf.gz
       bgzip varscan.cons
        """
}


process concatSnpsAndIndels {
    container = 'biocontainers/bcftools:v1.9-1-deb_cv1'

    input:
        tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi'), path('varscan.indels.vcf.gz'), path('varscan.indels.vcf.gz.tbi') from varscan_to_concat_qch

    output:
	tuple val(sampleName), path('varscan.concat.vcf') into varscan_concat_vcf_qch

    script:
	    """
        bcftools concat -a -o varscan.concat.vcf varscan.snps.vcf.gz varscan.indels.vcf.gz
	    """
}


process makeCombinedVarscanIndex {
    container = 'veupathdb/dnaseqanalysis'

    input:
	    tuple val(sampleName), path('varscan.concat.vcf') from varscan_concat_vcf_qch

    output:
	    tuple val(sampleName), path('varscan.concat.vcf.gz'), path('varscan.concat.vcf.gz.tbi') into varscan_concat_to_consensus_qch
        path('varscan.concat.vcf.gz') into (varscan_vcf_to_count_qch, varscan_vcf_to_merge_qch)
        path('varscan.concat.vcf.gz.tbi') into varscan_vcftbi_to_merge_qch

    script:
	    """
        bgzip varscan.concat.vcf
        tabix -fp vcf varscan.concat.vcf.gz
	    """
}


process mergeVcfs {
    container = 'biocontainers/bcftools:v1.9-1-deb_cv1'

    input:
        val 'vcfCount' from varscan_vcf_to_count_qch.collect().size()
        path('*.vcf.gz') from varscan_vcf_to_merge_qch.collect()
        path('*.vcf.gz.tbi') from varscan_vcftbi_to_merge_qch.collect()

    output:
	    path('result.vcf.gz') into merged_vcf_ch
    
    script:
        if (vcfCount > 1)
            """
            echo $vcfCount
            bcftools merge -o result.vcf.gz -O z *.vcf.gz
            """

        else
            """
            echo $vcfCount
            mv .vcf.gz result.vcf.gz
            """
}


process makeMergedVarscanIndex {
    container = 'veupathdb/dnaseqanalysis'

    publishDir "$params.outputDir", mode: "copy"

    input:
	    path('result.vcf.gz') from merged_vcf_ch

    output:
	    tuple path('result.vcf.gz'), path('result.vcf.gz.tbi')

    script:
	    """
        tabix -fp vcf result.vcf.gz
	    """
}


//TODO JB I have published this - please check
process bcftoolsConsensus {
    container = 'biocontainers/bcftools:v1.9-1-deb_cv1'

    publishDir "$params.outputDir", mode: "copy", saveAs: { filename -> "consensus.fa.gz" }

    input:
        tuple val(sampleName), path('varscan.concat.vcf.gz'), path('varscan.concat.vcf.gz.tbi') from varscan_concat_to_consensus_qch
        path 'genome.fa' from genome_fasta_path_ch
        path 'genome.fa.fai' from genome_index_path_ch

    output:
        path('cons.fa.gz')
    
    script:
	    """
        bcftools consensus -I -f genome.fa varscan.concat.vcf.gz > cons.fa
        gzip cons.fa
	    """
}


process genomecov {
    container = 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

    input:
        tuple val(sampleName), path('result_sorted_gatk.bam'), path('result_sorted_gatk.bai') from gatk_bam_to_genomecov_qch
        path 'genome.fa.fai' from genome_index_path_ch

    output:
	    tuple val(sampleName), path('coverage.bed') into bedgraph_qch

    script:
	    """
        bedtools genomecov -bg -ibam result_sorted_gatk.bam -g genome.fa.fai >coverage.bed
	    """
}


process bedGraphToBigWig {
    container = 'veupathdb/shortreadaligner'

    publishDir "$params.outputDir", mode: "copy", saveAs: { filename -> "${sampleName}.bw" }

    input:
        path 'genome.fa.fai' from genome_index_path_ch
        tuple val(sampleName), path('coverage.bed') from bedgraph_qch

    output:
	    path('coverage.bw')
    
    script:
	    """
        LC_COLLATE=C sort -k1,1 -k2,2n coverage.bed > sorted.coverage.bed
        bedGraphToBigWig sorted.coverage.bed genome.fa.fai coverage.bw
	    """
}


//CNV processes
process sortForCounting {
    container = 'veupathdb/shortreadaligner'

    input:
        tuple val(sampleName), path('result_sorted_gatk.bam'), path('result_sorted_gatk.bam.bai') from gatk_to_sortForCounting_qch

    output:
        tuple val(sampleName), path('result_sortByName.bam') into htseq_path_qch

    script:
        """
        samtools sort -n result_sorted_gatk.bam > result_sortByName.bam
        """
}


process htseqCount {
    container = 'biocontainers/htseq:v0.11.2-1-deb-py3_cv1'

    publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}.counts" }

    input:
        tuple val(sampleName), path('result_sortByName.bam') from htseq_path_qch
        path('gtfFile') from params.gtfFile

    output: 
        tuple val(sampleName), path('counts.txt') into count_to_tpm_qch

    script:
        """
        htseq-count -f bam -s no -t CDS -i gene_id -a 0 result_sortByName.bam gtfFile > counts.txt
        """
}


process calculateTPM {
    container = 'veupathdb/shortreadaligner'

    publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}.tpm" }

    input:
        tuple val(sampleName), path('counts.txt') from count_to_tpm_qch
        path('geneFootprintFile') from params.geneFootprintFile

    output:
        path('tpm.txt')

    script:
        """
        makeTpmFromHtseqCountsCNV.pl --geneFootprintFile geneFootprintFile --countFile counts.txt --tpmFile tpm.txt
        #NOTE downstream processing from here requires querying DBs and occasional reloading - leave in ReFlow
        """
}


process makeWindowFile {
    container = 'veupathdb/shortreadaligner'

    input:
        path('genome.fa.fai') from genome_index_path_ch
        val(winLen) from params.winLen

    output:
        tuple path('windows.bed'), path('genome.txt') into (window_to_bedtoolsWindowed_qch, window_to_snpDensity_qch, window_to_LOH_qch)

    script:
        """
        makeWindowedBed.pl --samtoolsIndex genome.fa.fai --winLen $winLen
        """
}


process bedtoolsWindowed {
    container = 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

    input:
        tuple path('windows.bed'), path('genome.txt') from window_to_bedtoolsWindowed_qch
        tuple val(sampleName), path('result_sorted_gatk.bam'), path('result_sorted_gatk.bam.bai') from gatk_to_bedtoolsWindowed_qch

    output:
        tuple val(sampleName), path('windowedCoverage.bed') into process_normaliseCoverage_qch

    script:
        """
        bedtools coverage -counts -sorted -g genome.txt -a windows.bed -b result_sorted_gatk.bam > windowedCoverage.bed
        """
}


process normaliseCoverage {
    container = 'veupathdb/shortreadaligner'

    publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}.bed" }

    input:
        tuple val(sampleName), path('windowedCoverage.bed'), path('summaryMetrics.txt') from process_normaliseCoverage_qch.join(picard_to_normaliseCoverage_qch)
   
    output:
        path('normalisedCoverage.bed')

    script:
        """
        # NOTE final processing requires querying the DB so can stay in ReFlow
        normaliseCoverageCNV.pl --bedFile windowedCoverage.bed --summaryMetrics summaryMetrics.txt
        """
} 


process makeSnpDensity {
    container= 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

    input:
        tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi'), path('varscan.indels.vcf.gz'), path('varscan.indels.vcf.gz.tbi') from varscan_to_snpDensity_qch
        tuple path('windows.bed'), path('genome.txt') from window_to_snpDensity_qch

    output:
        tuple val(sampleName), path('snpDensity.bed'), path('indelDensity.bed') into process_snpDensity_qch

    script: 
        """
        zcat varscan.snps.vcf.gz | bedtools coverage -a windows.bed -b stdin -sorted -g genome.txt -counts > snpDensity.bed
        zcat varscan.indels.vcf.gz | bedtools coverage -a windows.bed -b stdin -sorted -g genome.txt -counts > indelDensity.bed
        """
}


process makeDensityBigwigs {
    container = 'veupathdb/shortreadaligner'

    publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}_${filename}" }

    input:
        tuple val(sampleName), path('snpDensity.bed'), path('indelDensity.bed') from process_snpDensity_qch
        path('genome.fa.fai') from genome_index_path_ch

    output:
        tuple path('snpDensity.bw'), path('indelDensity.bw')

    script:
        """
        LC_COLLATE=C sort -k1,1 -k2,2n snpDensity.bed > sorted.snpDensity.bed
        LC_COLLATE=C sort -k1,1 -k2,2n indelDensity.bed > sorted.indelDensity.bed
        bedGraphToBigWig sorted.snpDensity.bed genome.fa.fai snpDensity.bw
        bedGraphToBigWig sorted.indelDensity.bed genome.fa.fai indelDensity.bw
        """
}

if(params.ploidy != 1) {
    process getHeterozygousSNPs {
        container = 'veupathdb/vcf_parser_cnv'
        input:
            tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi') from varscan_to_LOH_qch

        output:
            tuple val(sampleName), path('heterozygousSNPs.vcf') into process_makeLOH_qch

        script:
            """
            makeHeterozygosityPlot.py --vcfFile varscan.snps.vcf.gz 
            """
    }

    process makeHeterozygousDensityBed {
        container = 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'
        input:
            tuple path('windows.bed'), path('genome.txt') from window_to_LOH_qch
            tuple val(sampleName), path('heterozygousSNPs.vcf') from process_makeLOH_qch

        output:
            tuple val(sampleName), path('heterozygousDensity.bed') into process_heterozygousDensity_qch

        script:
            """
            bedtools coverage -a windows.bed -b heterozygousSNPs.vcf -sorted -g genome.txt -counts > heterozygousDensity.bed
            """
    }

    process makeHeterozygousDensityBigwig {
        container = 'veupathdb/shortreadaligner'

        publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}_LOH.bw" }

        input:
            tuple val(sampleName), path('heterozygousDensity.bed') from process_heterozygousDensity_qch
            path('genome.fa.fai') from genome_index_path_ch

        output:
            path('heterozygousDensity.bw') 

        script:
            """
            LC_COLLATE=C sort -k1,1 -k2,2n heterozygousDensity.bed > sorted.heterozygousDensity.bed
            bedGraphToBigWig sorted.heterozygousDensity.bed genome.fa.fai heterozygousDensity.bw
            """
    }
}
