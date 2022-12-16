# <p align=center>dnaseqAnalysis Nextflow Workflows [^1]</p>

[^1]: Both workflows are still in development   

***<p align=center>processSingleExperiment</p>***  
```mermaid
flowchart TD
    p0((Channel.fromPath))
    p1([map])
    p2(( ))
    p3[processSingleExperiment:hisat2Index]
    p4[processSingleExperiment:fastqc]
    p5([join])
    p6[processSingleExperiment:fastqc_check]
    p7([join])
    p8(( ))
    p9[processSingleExperiment:trimmomatic]
    p10([join])
    p11([join])
    p12[processSingleExperiment:hisat2]
    p13([first])
    p14(( ))
    p15[processSingleExperiment:reorderFasta]
    p16[processSingleExperiment:subsample]
    p17(( ))
    p18[processSingleExperiment:picard]
    p19(( ))
    p20[processSingleExperiment:gatk]
    p21[processSingleExperiment:mpileup]
    p22([join])
    p23(( ))
    p24[processSingleExperiment:varscan]
    p25(( ))
    p26[processSingleExperiment:concatSnpsAndIndels]
    p27[processSingleExperiment:makeCombinedVarscanIndex]
    p28[processSingleExperiment:filterIndels]
    p29[processSingleExperiment:makeIndelTSV]
    p30(( ))
    p31([count])
    p32([map])
    p33([groupTuple])
    p34[processSingleExperiment:mergeVcfs]
    p35[processSingleExperiment:makeMergedVarscanIndex]
    p36(( ))
    p37[processSingleExperiment:bcftoolsConsensus]
    p38[processSingleExperiment:addSampleToDefline]
    p39(( ))
    p40[processSingleExperiment:genomecov]
    p41[processSingleExperiment:bedGraphToBigWig]
    p42(( ))
    p43[processSingleExperiment:sortForCounting]
    p44(( ))
    p45[processSingleExperiment:htseqCount]
    p46(( ))
    p47[processSingleExperiment:calculateTPM]
    p48(( ))
    p49(( ))
    p50[processSingleExperiment:makeWindowFile]
    p51[processSingleExperiment:bedtoolsWindowed]
    p52([join])
    p53[processSingleExperiment:normaliseCoverage]
    p54(( ))
    p55[processSingleExperiment:makeSnpDensity]
    p56[processSingleExperiment:makeDensityBigwigs]
    p57(( ))
    p0 --> p1
    p1 -->|samples_qch| p4
    p2 -->|genome.fa| p3
    p3 --> p12
    p3 --> p12
    p4 --> p5
    p1 -->|samples_qch| p5
    p5 --> p6
    p6 --> p7
    p1 -->|samples_qch| p7
    p7 --> p9
    p8 -->|adaptorsFile| p9
    p9 --> p11
    p1 -->|samples_qch| p10
    p6 --> p10
    p10 --> p11
    p11 --> p12
    p12 --> p13
    p13 --> p15
    p14 -->|genome.fa| p15
    p15 --> p18
    p12 --> p16
    p16 --> p18
    p17 -->|picardJar| p18
    p18 --> p20
    p18 --> p52
    p19 -->|gatkJar| p20
    p15 --> p20
    p20 --> p21
    p15 --> p21
    p21 --> p22
    p20 --> p22
    p22 --> p24
    p23 -->|varscanJar| p24
    p15 --> p24
    p24 --> p26
    p24 --> p25
    p26 --> p27
    p27 --> p28
    p28 --> p29
    p29 --> p30
    p27 --> p31
    p31 --> p34
    p27 --> p32
    p32 --> p33
    p33 --> p34
    p34 --> p35
    p35 --> p36
    p27 --> p37
    p15 --> p37
    p37 --> p38
    p38 --> p39
    p20 --> p40
    p15 --> p40
    p40 --> p41
    p15 --> p41
    p41 --> p42
    p20 --> p43
    p43 --> p45
    p44 -->|gtfFile| p45
    p45 --> p47
    p46 -->|geneFootprintFile| p47
    p47 --> p48
    p15 --> p50
    p49 -->|winLen| p50
    p50 --> p51
    p20 --> p51
    p51 --> p52
    p52 --> p53
    p53 --> p54
    p24 --> p55
    p50 --> p55
    p55 --> p56
    p15 --> p56
    p56 --> p57
```
This workflow will run on a per organism basis with multiple strains. Anyone can run this workflow, as it does not require a gus environment. The output from this process will either be sent to webservices, uploaded to various databases, and/or used in the mergeExperiments process.  
  

***<p align=center>loadCNV</p>***  
```mermaid
flowchart TD
    p0((Channel.fromPath))
    p1([map])
    p2(( ))
    p3(( ))
    p4(( ))
    p5(( ))
    p6[loadCNV:calculatePloidyAndGeneCNV]
    p7(( ))
    p8[loadCNV:writePloidyConfigFile]
    p9[loadCNV:writeCNVConfigFile]
    p10[loadCNV:loadPloidy]
    p11[loadCNV:loadGeneCNV]
    p0 --> p1
    p1 -->|tpm_qch| p6
    p2 -->|footprints| p6
    p3 -->|gusConfig| p6
    p4 -->|ploidy| p6
    p5 -->|taxonId| p6
    p6 --> p8
    p6 --> p9
    p6 --> p7
    p8 --> p10
    p8 -->|ploidyFile| p10
    p9 --> p11
    p9 -->|geneCNVFile| p11
```

***<p align=center>loadSingleExperiment</p>***  
```mermaid
flowchart TD
    p0((Channel.fromPath))
    p1((Channel.fromPath))
    p2((Channel.fromPath))
    p3(( ))
    p4(( ))
    p5[loadSingleExperiment:loadIndels]
    p6([collectFile])
    p7(( ))
    p8([collectFile])
    p9(( ))
    p0 -->|indels_qch| p5
    p1 -->|bam_qch| p8
    p2 -->|bw_qch| p6
    p3 -->|extDbRlsSpec| p5
    p4 -->|genomeExtDbRlsSpec| p5
    p6 --> p7
    p8 --> p9
```
This work flow will run after processSingleExperiment. This workflow requires a gus environment to run. It will take the strain specific vcfs and consensus sequences output from the processSingleExperiment workflow. The strain specific vcfs will be merged together to create a merged vcf. This will be sent to webservices, along with being sent to snpEff to generate an annotated vcf file. The consensus sequences will be combined and sent to web services. These masked consensus sequences, along with various information queried from our databases will be used to generate a transcript fasta file that is indel and coverage aware. This will be used in downstream processes that still need to be generated.

***<p align=center>mergeExperiments</p>***  

```mermaid
flowchart TD
    p0((Channel.fromPath))
    p1((Channel.fromPath))
    p2([collectFile])
    p3[mergeExperiments:checkUniqueIds]
    p4(( ))
    p5([collect])
    p6[mergeExperiments:mergeVcfs]
    p7(( ))
    p8[mergeExperiments:makeSnpFile]
    p9(( ))
    p10(( ))
    p11(( ))
    p12(( ))
    p13(( ))
    p14(( ))
    p15(( ))
    p16(( ))
    p17(( ))
    p18[mergeExperiments:processSeqVars]
    p19(( ))
    p20(( ))
    p21(( ))
    p22(( ))
    p23(( ))
    p24(( ))
    p25[mergeExperiments:snpEff]
    p26(( ))
    p0 -->|fastas_qch| p2
    p1 -->|vcfs_qch| p5
    p2 -->|combinedFastagz| p3
    p3 -->|-| p4
    p5 -->|allvcfs| p6
    p6 --> p7
    p6 --> p8
    p8 --> p18
    p9 -->|gusConfig.txt| p18
    p10 -->|cache.txt| p18
    p11 -->|undoneStrains.txt| p18
    p12 -->|transcript_extdb_spec| p18
    p13 -->|organism_abbrev| p18
    p14 -->|reference_strain| p18
    p15 -->|extdb_spec| p18
    p16 -->|varscan_directory| p18
    p17 -->|genome.fa| p18
    p2 -->|combinedFastagz| p18
    p18 --> p22
    p18 --> p21
    p18 --> p20
    p18 --> p19
    p6 -->|merged.vcf| p25
    p23 -->|genes.gtf.gz| p25
    p24 -->|sequences.fa.gz| p25
    p25 --> p26

```

**<p align=center>Explanation of Config File Parameters</p>**
---

Workflows: processSingleExperiment=ps; loadSingleExperiment=ls; loadCNV=lc; mergeExperiments=me;

| ps | ls | lc | me | Parameter | Value | Description |
| -- | -- | -- | -- | --------- | ----- | ----------- |
| X |   | X | X | outputDir | string path | Where you would like the output files to be stored |
| X | X | X |   | inputDir | string path | Path to the directory containing the strain specific fastqs and bam files (each strain has their own directory located inside the inputDir) |
| X |   |   |   | fromBAM | boolean | If true, samples will be retrieved from the strain specific bam files |
| X |   |   |   | hisat2Threads | int | Specifies NTHREADS parallel search threads count, to be used as argument -p when calling hisat2 |
| X |   |   |   | isPaired | boolean | (Assuming fromBam is false) Specifies if the samples are being retrieved from paired or unpaired fastq files |
| X |   |   |   | minCoverage | int | Sets the minimum coverage value, used in varscan (mpileup2snp, mpileup2indel, mpileup2cns) as well as the minimum coverage value used to determine masking when creating the consensus genome per strain |
| X |   |   | X | genomeFastaFile | string path | Path to the genome fasta file to be used as reference. Also used in used in snpEff database creation |
| X |   |   | X | gtfFile | string path | Path to the gtf file. Used in snpEff database creation (me) |
| X |   |   |   | winLen | int | Specifies the window length argument used in bin/makeWindowedBed.pl |
| X |   |   | X | ploidy | int | Ploidy Level |
| X |   |   |   | hisat2Index | string path | (Assuming createIndex is false) Location of the hisat2Index file  |
| X |   |   |   | createIndex | boolean | If true, will create the hisat2Index file | 
| X |   |   |   | trimmomaticAdaptorsFile | string path | Location of the trimmomatic adaptors file |
| X |   |   |   | varscanPValue | int | Sets the --p-value argument used in varscan mpileup2snp, mpileup2indel, and mpileup2cns |
| X |   |   |   | varscanMinVarFreqSnp | int | Sets the --min-var-freq argument used in varscan mpileup2snp |
| X |   |   |   | varscanMinVarFreqCons | int | Sets the --min-var-freq argument used in varscan mpileup2indel and mpileup2cns |
| X |   |   |   | maxNumberOfReads | int | Used in subSample process to limit total number of reads | 
| X |   | X |   | footprintFile | path | Path to gene footprints file |
|   |   |   | X | fastaDir | string path | Path to directory that contains the consensus fasta files output from processSingleExperiment |
|   |   |   | X | vcfDir | string path | Path to directory that contains the strain specific vcf files output from processSingleExperiment |
|   |   |   | X | makepositionarraycoding | string path | Path to makePotionArrayCoding.pl Leave as 'bin/makePositionArrayCoding', this parameter will most likely be removed in later updates |
|   |   | X | X | gusConfig | path | Path to gus.config file |
|   |   |   | X | cacheFile | path | Path to cache file |
|   |   |   | X | undoneStrains | path | Path to undoneStrains file |
|   |   |   | X | organism_abbrev | string | Organism Abbreviation Ex: 'lmajFriedlin' |
|   |   |   | X | reference_strain | string | Reference Strain Ex: 'Friedlin' |
|   |   |   | X | varscan_directory | path | Path to varscan coverage directory |
|   |   | X |   | taxonId | string | Taxon ID Ex: "1577702" |
|   | X |   | X | extDbRlsSpec | string | External database release spec. Ex: "lmajFriedlin_NGS_SNPsAndVariations|do_not_care" |
|   | X |   | X | genomeExtDbRlsSpec | string | Genome external database spec. Ex: "lmajFriedlin_primary_genome_RSRC|2016-05-28" |
|   | X |   |   | webServicesDir | path | Path to web services directory |
