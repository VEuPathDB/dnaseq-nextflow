# dnaSeqAnalysis Nextflow

**Explanation of Workflows:**


**Explanation of Config File Parameters:**

| Workflow | Parameter | Value | Description |
| -------- | --------- | ----- | ----------- |
| processSingleExperiment & mergeExperiments | workflow | string | Either 'processSingleExperiment' or 'mergeExperiments', depending on what workflow you would like to run |
| processSingleExperiment & mergeExperiments | outputDir | string path | Where you would like the output files to be stored |
| processSingleExperiment | inputDir | string path | Path to the directory containing the strain specific fastqs and bam files (each strain has their own directory located inside the inputDir) |
| processSingleExperiment | fromBAM | boolean | If true, samples will be retrieved from the strain specific bam files |
| processSingleExperiment | hisat2Threads | int | Specifies NTHREADS parallel search threads count, to be used as argument -p when calling hisat2 |
| processSingleExperiment | isPaired | boolean | (Assuming fromBam is false) Specifies if the samples are being retrieved from paired or unpaired fastq files |
| processSingleExperiment | minCoverage | int | Sets the minimum coverage value, used in varscan (mpileup2snp, mpileup2indel, mpileup2cns) as well as the minimum coverage value used to determine masking when creating the consensus genome per strain |
| processSingleExperiment | genomeFastaFile | string path | Path to the genome fasta file to be used as reference |
| processSingleExperiment | gtfFile | string path | Path to the gtf file  |
| processSingleExperiment | geneFootprintFile | string path | Path to the footprint file |
| processSingleExperiment | winLen | int | Specifies the window length argument used in bin/makeWindowedBed.pl |
| processSingleExperiment | ploidy | int | Ploidy Level |
| processSingleExperiment | hisat2Index | string path | (Assuming createIndex is false) Location of the hisat2Index file  |
| processSingleExperiment | createIndex | boolean | If true, will create the hisat2Index file | 
| processSingleExperiment | trimmomaticAdaptorsFile | string path | Location of the trimmomatic adaptors file |
| processSingleExperiment | varscanPValue | int | Sets the --p-value argument used in varscan mpileup2snp, mpileup2indel, and mpileup2cns |
| processSingleExperiment | varscanMinVarFreqSnp | int | Sets the --min-var-freq argument used in varscan mpileup2snp |
| processSingleExperiment | varscanMinVarFreqCons | int | Sets the --min-var-freq argument used in varscan mpileup2indel and mpileup2cns |
| processSingleExperiment | maxNumberOfReads | int | Used in subSample process to limit total number of reads | 
| processSingleExperiment | gatkJar | string path | Path to gatk Jar. If defined as 'NA', jar path will be defined to work inside of the container |
| processSingleExperiment | picardJar | string path | Path to picard Jar. If defined as 'NA', jar path will be defined to work inside of the container  |
| processSingleExperiment | varscanJar | string path | Path to varscan Jar. If defined as 'NA', jar path will be defined to work inside of the container  |
| mergeExperiments | fastaDir | string path | Path to directory that contains the consensus fasta files output from processSingleExperiment |
| mergeExperiments | vcfDir | string path | Path to directory that contains the strain specific vcf files output from processSingleExperiment |
| mergeExperiments | makepositionarraycoding | string path | Path to makePotionArrayCoding.pl Leave as 'bin/makePositionArrayCoding', this parameter will most likely be removed in later updates |
| mergeExperiments | databaseFileType | string | Either 'gtf' or 'gff', depending on the annotation file type to be used in snpEff database creation |
| mergeExperiments | databaseFile | string path | Path to the gtf or gff file to be used in snpEff database creation |
| mergeExperiments | sequenceFile | string path | Path to the reference genome fasta to be used in snpEff database creation |
