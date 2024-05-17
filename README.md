# <p align=center>dnaseqAnalysis Nextflow Workflows [^1]</p>

[^1]: Both workflows are still in development   

***<p align=center>processSingleExperiment</p>***  
```mermaid
flowchart TD
    p0((Channel.fromFilePairs))
    p1(( ))
    p2(( ))
    p3(( ))
    p4[processSingleExperiment:ps:hisat2Index]
    p5(( ))
    p6[processSingleExperiment:ps:fastqc]
    p7([join])
    p8(( ))
    p9[processSingleExperiment:ps:fastqc_check]
    p10([join])
    p11(( ))
    p12(( ))
    p13[processSingleExperiment:ps:trimmomatic]
    p14([join])
    p15([join])
    p16(( ))
    p17(( ))
    p18[processSingleExperiment:ps:hisat2]
    p19([first])
    p20(( ))
    p21[processSingleExperiment:ps:reorderFasta]
    p22[processSingleExperiment:ps:subsample]
    p23[processSingleExperiment:ps:picard]
    p24[processSingleExperiment:ps:gatk]
    p25[processSingleExperiment:ps:mpileup]
    p26([join])
    p27[processSingleExperiment:ps:varscan]
    p28(( ))
    p29[processSingleExperiment:ps:concatSnpsAndIndels]
    p30[processSingleExperiment:ps:makeCombinedVarscanIndex]
    p31[processSingleExperiment:ps:filterIndels]
    p32[processSingleExperiment:ps:makeIndelTSV]
    p33(( ))
    p34([count])
    p35([map])
    p36([groupTuple])
    p37[processSingleExperiment:ps:mergeVcfs]
    p38[processSingleExperiment:ps:makeMergedVarscanIndex]
    p39(( ))
    p40[processSingleExperiment:ps:bcftoolsConsensus]
    p41[processSingleExperiment:ps:addSampleToDefline]
    p42(( ))
    p43[processSingleExperiment:ps:genomecov]
    p44[processSingleExperiment:ps:bedGraphToBigWig]
    p45(( ))
    p46[processSingleExperiment:ps:sortForCounting]
    p47(( ))
    p48[processSingleExperiment:ps:htseqCount]
    p49(( ))
    p50[processSingleExperiment:ps:calculateTPM]
    p51(( ))
    p52(( ))
    p53(( ))
    p54(( ))
    p55(( ))
    p56[processSingleExperiment:ps:calculatePloidyAndGeneCNV]
    p57(( ))
    p58(( ))
    p59(( ))
    p60(( ))
    p61[processSingleExperiment:ps:makeWindowFile]
    p62[processSingleExperiment:ps:bedtoolsWindowed]
    p63([join])
    p64[processSingleExperiment:ps:normaliseCoverage]
    p65(( ))
    p66[processSingleExperiment:ps:makeSnpDensity]
    p67[processSingleExperiment:ps:makeDensityBigwigs]
    p68(( ))
    p69[processSingleExperiment:ps:getHeterozygousSNPs]
    p70[processSingleExperiment:ps:makeHeterozygousDensityBed]
    p71[processSingleExperiment:ps:makeHeterozygousDensityBigwig]
    p72(( ))
    p0 -->|samples_qch| p6
    p1 -->|genomeFasta| p4
    p2 -->|fromBam| p4
    p3 -->|createIndex| p4
    p4 --> p18
    p4 --> p18
    p5 -->|fromBam| p6
    p6 --> p7
    p0 -->|samples_qch| p7
    p7 --> p9
    p8 -->|fromBam| p9
    p9 --> p10
    p0 -->|samples_qch| p10
    p10 --> p13
    p11 -->|fromBam| p13
    p12 -->|isPaired| p13
    p13 --> p15
    p9 --> p14
    p0 -->|samples_qch| p14
    p14 --> p15
    p15 --> p18
    p16 -->|fromBam| p18
    p17 -->|isPaired| p18
    p18 --> p19
    p19 --> p21
    p20 -->|genomeFasta| p21
    p21 --> p23
    p18 --> p22
    p22 --> p23
    p23 --> p24
    p23 --> p63
    p21 --> p24
    p24 --> p25
    p21 --> p25
    p25 --> p26
    p24 --> p26
    p26 --> p27
    p21 --> p27
    p27 --> p29
    p27 --> p28
    p29 --> p30
    p30 --> p31
    p31 --> p32
    p32 --> p33
    p30 --> p34
    p34 --> p37
    p30 --> p35
    p35 --> p36
    p36 --> p37
    p37 --> p38
    p38 --> p39
    p30 --> p40
    p21 --> p40
    p40 --> p41
    p41 --> p42
    p24 --> p43
    p21 --> p43
    p43 --> p44
    p21 --> p44
    p44 --> p45
    p24 --> p46
    p46 --> p48
    p47 -->|gtfFile| p48
    p48 --> p50
    p49 -->|geneFootprintFile| p50
    p50 --> p56
    p51 -->|footprints| p56
    p52 -->|ploidy| p56
    p53 -->|taxonId| p56
    p54 -->|geneSourceIdOrtholog| p56
    p55 -->|chrsForCalc| p56
    p56 --> p59
    p56 --> p58
    p56 --> p57
    p21 --> p61
    p60 -->|winLen| p61
    p61 --> p62
    p24 --> p62
    p62 --> p63
    p63 --> p64
    p64 --> p65
    p27 --> p66
    p61 --> p66
    p66 --> p67
    p21 --> p67
    p67 --> p68
    p27 --> p69
    p69 --> p70
    p61 --> p70
    p70 --> p71
    p21 --> p71
    p71 --> p72
```
This workflow will run on a per organism basis with multiple strains. Anyone can run this workflow, as it does not require a gus environment. The output from this process will either be sent to webservices, uploaded to various databases, and/or used in the mergeExperiments process.  
  
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

| ps | me | Parameter | Value | Description |
| -- | -- | -- | -- | --------- |
| X  | X | outputDir | string path | Where you would like the output files to be stored |
| X  |   | inputDir | string path | Path to the directory containing the strain specific fastqs and bam files (each strain has their own directory located inside the inputDir) |
| X  |   | fromBAM | boolean | If true, samples will be retrieved from the strain specific bam files |
| X  |   | hisat2Threads | int | Specifies NTHREADS parallel search threads count, to be used as argument -p when calling hisat2 |
| X  |   | isPaired | boolean | (Assuming fromBam is false) Specifies if the samples are being retrieved from paired or unpaired fastq files |
| X  |   | minCoverage | int | Sets the minimum coverage value, used in varscan (mpileup2snp, mpileup2indel, mpileup2cns) as well as the minimum coverage value used to determine masking when creating the consensus genome per strain |
| X | X | genomeFastaFile | string path | Path to the genome fasta file to be used as reference. Also used in used in snpEff database creation |
| X | X | gtfFile | string path | Path to the gtf file. Used in snpEff database creation (me) |
| X |   | winLen | int | Specifies the window length argument used in bin/makeWindowedBed.pl |
| X | X | ploidy | int | Ploidy Level |
| X |   | hisat2Index | string path | (Assuming createIndex is false) Location of the hisat2Index file  |
| X |   | createIndex | boolean | If true, will create the hisat2Index file | 
| X |   | trimmomaticAdaptorsFile | string path | Location of the trimmomatic adaptors file |
| X |   | varscanPValue | int | Sets the --p-value argument used in varscan mpileup2snp, mpileup2indel, and mpileup2cns |
| X |   | varscanMinVarFreqSnp | int | Sets the --min-var-freq argument used in varscan mpileup2snp |
| X |   | varscanMinVarFreqCons | int | Sets the --min-var-freq argument used in varscan mpileup2indel and mpileup2cns |
| X |   | maxNumberOfReads | int | Used in subSample process to limit total number of reads | 
| X |   | footprintFile | path | Path to gene footprints file |
|   | X | fastaDir | string path | Path to directory that contains the consensus fasta files output from processSingleExperiment |
|   | X | vcfDir | string path | Path to directory that contains the strain specific vcf files output from processSingleExperiment |
|   | X | makepositionarraycoding | string path | Path to makePotionArrayCoding.pl Leave as 'bin/makePositionArrayCoding', this parameter will most likely be removed in later updates |
|   | X | gusConfig | path | Path to gus.config file |
|   | X | cacheFile | path | Path to cache file |
|   | X | undoneStrains | path | Path to undoneStrains file |
|   | X | organism_abbrev | string | Organism Abbreviation Ex: 'lmajFriedlin' |
|   | X | reference_strain | string | Reference Strain Ex: 'Friedlin' |
|   | X | varscan_directory | path | Path to varscan coverage directory |
|   |   | taxonId | string | Taxon ID Ex: "1577702" |
|   | X | extDbRlsSpec | string | External database release spec. Ex: "lmajFriedlin_NGS_SNPsAndVariations|do_not_care" |
|   | X | genomeExtDbRlsSpec | string | Genome external database spec. Ex: "lmajFriedlin_primary_genome_RSRC|2016-05-28" |
|   |   | webServicesDir | path | Path to web services directory |
