if(params.fromBAM) {
    samples_ch = Channel.fromPath([params.inputDir + '/**/*.bam']).map { file -> tuple(file.baseName, [file]) }
}
else if(params.isPaired) {
   samples_ch = Channel.fromFilePairs([params.inputDir + '/**/*_{1,2}.fastq', params.inputDir + '/**/*_{1,2}.fastq.gz'])
}
else {
    samples_ch = Channel.fromPath([params.inputDir + '/**/*.fastq', params.inputDir + '/**/*.fastq.gz']).map { file -> tuple(file.baseName, [file]) }
}

gatkJar = file(params.gatkJar);
picardJar = file(params.picardJar);
varscanJar = file(params.varscanJar);

process hisat2Index {
     input:
     	path 'genome.fa' from params.genomeFastaFile
    
    output:
	path 'genomeIndex*.ht2' into hisat2_index_file_ch
        val 'genomeIndex' into hisat2_index_ch
        path 'genome.fa.fai' into genome_index_ch
        path 'genome.fa' into genome_fasta_ch

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
  input:
	tuple val(sampleName), path(sampleFile) from samples_ch

  output:
	tuple val(sampleName), path(sampleFile) into fastqc_samples_ch
        path('fastqc_output', type:'dir') into fastqc_dir_ch

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
  input:
	tuple val(sampleName), path(sampleFile) from fastqc_samples_ch
        path 'fastqc_output' from fastqc_dir_ch

  output:
	tuple val(sampleName), path(sampleFile) into fastqc_check_samples_ch
        path 'mateAEncoding' into fastqc_check_encoding_ch

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
  input:
	tuple val(sampleName), path(sampleFile) from fastqc_check_samples_ch
        path('mateAEncoding') from fastqc_check_encoding_ch
        path adaptorsFile from params.trimmomaticAdaptorsFile

  output:
	tuple val(sampleName), path(sampleFile) into trimmomatic_samples_ch
        path 'mateAEncoding' into trimmomatic_encoding_ch
        path 'sample_1P' into trimmomatic_1p_ch
        path 'sample_2P' into trimmomatic_2p_ch

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
        touch sample_2p
        mateAEncoding=\$(<mateAEncoding)
        java org.usadellab.trimmomatic.TrimmomaticSE -trimlog trimLog.txt -$mateAEncoding $sampleFile sample ILLUMINACLIP:$adaptorsFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
        """
}


process hisat2 {
  input:
	tuple val(sampleName), path(sampleFile) from trimmomatic_samples_ch
        path 'mateAEncoding' from trimmomatic_encoding_ch
        path 'sample_1p' from trimmomatic_1p_ch
        path 'sample_2p' from trimmomatic_2p_ch
        val hisat2_index from hisat2_index_ch
        path 'genomeIndex.*.ht2' from hisat2_index_file_ch

    
    output:
	val(sampleName) into hisat2_samples_ch
        path('result_sorted.bam') into hisat2_result_bam_ch
    
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

process picard {
     input:
        path 'genome.fa' from genome_fasta_ch
        path 'genome.fa.fai' from genome_index_ch
        path picardJar from gatkJar
	val(sampleName) from hisat2_samples_ch
        path('result_sorted.bam') from hisat2_result_bam_ch


     output:
        path 'genome.dict' into picard_dict_ch
	val(sampleName) into picard_samples_ch
        path('picard.bam') into picard_result_bam_ch
        path('picard.bai') into picard_result_bai_ch


    script:
        def jarPath = picardJar.name == 'NA' ? "/usr/picard/picard.jar" : picardJar.name
        """
        # TODO Mapping stats and Downsample
        java -jar $jarPath AddOrReplaceReadGroups I=result_sorted.bam O=picard.bam RGID=$sampleName RGSM=$sampleName RGLB=NA RGPL=NA RGPU=NA
        java -jar $jarPath BuildBamIndex I=picard.bam
        java -jar $jarPath CreateSequenceDictionary R=genome.fa UR=genome.fa
        """


}

process gatk {
    input:
        path gatkJar from gatkJar
        path 'genome.fa' from genome_fasta_ch
        path 'genome.fa.fai' from genome_index_ch
        path 'genome.dict' from picard_dict_ch
	val(sampleName) from picard_samples_ch
        path('picard.bam') from picard_result_bam_ch
        path('picard.bai') from picard_result_bai_ch

    output:
	val(sampleName) into (gatk_samples_ch, genome_coverage_samples_ch)
        path('result_sorted_gatk.bam') into (gatk_result_bam_ch, genome_coverage_bam_ch)
        path('result_sorted_gatk.bai') into (gatk_result_bai_ch, genome_coverage_bai_ch)

    script:
        def jarPath = gatkJar.name == 'NA' ? "/usr/GenomeAnalysisTK.jar" : gatkJar.name

        """
        java -jar $jarPath -I picard.bam -R genome.fa -T RealignerTargetCreator -o forIndelRealigner.intervals 2>realaligner.err
        java -jar $jarPath -I picard.bam -R genome.fa -T IndelRealigner -targetIntervals forIndelRealigner.intervals -o result_sorted_gatk.bam 2>indelRealigner.err
        """
}


process mpileup {
    input:
	val(sampleName) from gatk_samples_ch
        path ('result_sorted_gatk.bam') from gatk_result_bam_ch
        path('result_sorted_gatk.bam.bai') from gatk_result_bai_ch
        path 'genome.fa' from genome_fasta_ch
        path 'genome.fa.fai' from genome_index_ch

    
    output:
	val(sampleName) into mpileup_samples_ch
        path ('result_sorted_gatk.bam') into mpileup_result_bam_ch
        path('result_sorted_gatk.bam.bai') into mpileup_result_bai_ch
        path('result.pileup') into mpileup_result_ch
    
    script:
	"""
        samtools mpileup -A -f genome.fa -B result_sorted_gatk.bam > result.pileup 2>pileup.err
        """
}


process varscan {
    input:
	val(sampleName) from mpileup_samples_ch
        path ('result_sorted_gatk.bam') from mpileup_result_bam_ch
        path('result_sorted_gatk.bam.bai') from mpileup_result_bai_ch
        path('result.pileup') from mpileup_result_ch

    output:
	val(sampleName) into varscan_samples_ch
        path('varscan.snps.vcf.gz') into varscan_snp_ch
        path('varscan.snps.vcf.gz.tbi') into varscan_snp_indx_ch
    
    
    script:
        def jarPath = varscanJar.name == 'NA' ? "/usr/local/VarScan.jar" : varscanJar.name
	"""
       java -jar $jarPath mpileup2snp result.pileup --output-vcf 1 --p-value $params.varscanPValue --min-coverage $params.varscanMinCoverage --min-var-freq $params.varscanMinVarFreqSnp >varscan.snps.vcf  2>varscan_snps.err
       java -jar $jarPath mpileup2indel result.pileup --output-vcf 1 --p-value $params.varscanPValue --min-coverage $params.varscanMinCoverage --min-var-freq $params.varscanMinVarFreqCons >varscan.indels.vcf  2> varscan_indels.err
       java -jar $jarPath mpileup2cns result.pileup --p-value $params.varscanPValue --min-coverage $params.varscanMinCoverage --min-var-freq $params.varscanMinVarFreqCons > varscan.cons 2>varscan_cons.err
       bgzip varscan.snps.vcf
       tabix -fp vcf varscan.snps.vcf.gz

        """
	
    
}


process bcftools {
    input:
	val(sampleName) from varscan_samples_ch
        path('varscan.snps.vcf.gz') from varscan_snp_ch
        path('varscan.snps.vcf.gz.tbi') from varscan_snp_indx_ch
        path 'genome.fa' from genome_fasta_ch
        path 'genome.fa.fai' from genome_index_ch

    
    script:
	"""
        bcftools consensus -I -f genome.fa varscan.snps.vcf.gz > cons.fa
        gzip cons.fa
	"""
}

process bedtools {
    input:
	val(sampleName) from genome_coverage_samples_ch
        path('result_sorted_gatk.bam') from genome_coverage_bam_ch
        path('result_sorted_gatk.bai') from genome_coverage_bai_ch
        path 'genome.fa.fai' from genome_index_ch

    output:
	path('coverage.bed') into bedgraph_ch
	val(sampleName) into bedgraph_samples_ch

    script:
	"""
        bedtools genomecov -bg -ibam result_sorted_gatk.bam -g genome.fa.fai >coverage.bed
	"""
}




process bedGraphToBigWig {
    input:
	val(sampleName) from bedgraph_samples_ch
        path 'genome.fa.fai' from genome_index_ch
	path('coverage.bed') from bedgraph_ch

    script:
	"""
        bedGraphToBigWig coverage.bed genome.fa.fai coverage.bw
	"""
}
