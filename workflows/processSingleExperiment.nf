#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Preprocessing
include { fastqc } from '../modules/preprocessing.nf'
include { fastqc_check } from '../modules/preprocessing.nf'
include { trimmomatic } from '../modules/preprocessing.nf'

// Alignment
include { bwaIndex } from '../modules/alignment.nf'
include { bwaMem } from '../modules/alignment.nf'
include { reorderFasta } from '../modules/alignment.nf'
include { picard } from '../modules/alignment.nf'
include { gatk } from '../modules/alignment.nf'
include { rawReadCount } from '../modules/alignment.nf'
include { samtoolsStats } from '../modules/alignment.nf'
include { bedtoolsGenomecovStats } from '../modules/alignment.nf'
include { makeAlignmentStats } from '../modules/alignment.nf'

// SNP
include { runFreebayes } from '../modules/snp.nf'
include { filterAndSplitVcf } from '../modules/snp.nf'
include { makeIndelTSV } from '../modules/snp.nf'
include { makeCoverageBed } from '../modules/snp.nf'
include { makeConsensusFromVcfAndBed } from '../modules/snp.nf'

// CNV
include { genomecov } from '../modules/cnv.nf'
include { bedGraphToBigWig } from '../modules/cnv.nf'
include { normaliseCoverageToBigWig } from '../modules/cnv.nf'

include { sortForCounting } from '../modules/cnv.nf'
include { htseqCount } from '../modules/cnv.nf'
include { calculateTPM } from '../modules/cnv.nf'
include { makeWindowFile } from '../modules/cnv.nf'
include { bedtoolsWindowed } from '../modules/cnv.nf'
include { normaliseCoverage } from '../modules/cnv.nf'
include { makeSnpDensity } from '../modules/cnv.nf'
include { makeDensityBigwigs } from '../modules/cnv.nf'
include { convertFreebayesToVarscanFormat } from '../modules/cnv.nf'
include { getHeterozygousSNPs } from '../modules/cnv.nf'
include { makeHeterozygousDensityBed } from '../modules/cnv.nf'
include { makeHeterozygousDensityBigwig } from '../modules/cnv.nf'
include { calculatePloidyAndGeneCNV } from '../modules/cnv.nf'

workflow ps {

  take:

    samples_qch

  main:

    genome_fasta_file = file(params.genomeFastaFile)

    // Contig-only assemblies can leave chrsForCalcFile empty (no chromosomes to
    // normalise against), which would make ploidy/gene-CNV estimation and the
    // normalised-coverage track fail or produce bogus output. Skip just that
    // functionality in that case; SNP calling and density tracks are unaffected.
    chrsForCalcFile = file(params.chrsForCalcFile)
    runCnv = chrsForCalcFile.exists() && chrsForCalcFile.text.trim().length() > 0

    bwaIndexResults = bwaIndex(genome_fasta_file)

    // Extract is_paired and files channels separately
    samples_qch
      .map { sample_id, files, is_paired -> tuple(sample_id, files) }
      .set { files_only_qch }

    samples_qch
      .map { sample_id, files, is_paired -> tuple(sample_id, is_paired) }
      .set { paired_info_qch }

    fastqcResults = fastqc(files_only_qch)

    fastqc_checkResults = fastqc_check(files_only_qch.join(fastqcResults))

    trimmomaticResults = trimmomatic(
      files_only_qch.join(fastqc_checkResults).join(paired_info_qch)
    )

    bwaMemResults = bwaMem(
      files_only_qch.join(fastqc_checkResults).join(trimmomaticResults).join(paired_info_qch),
      bwaIndexResults.index_files,
      bwaIndexResults.genome_fasta
    )

    reorderFastaResults = reorderFasta(bwaMemResults.first(), genome_fasta_file)


    picardResults = picard(reorderFastaResults, bwaMemResults)

    gatkResults = gatk(reorderFastaResults, picardResults.bam_and_dict )

    genomeRepeatsBed = params.genomeRepeatsBedFile ? file(params.genomeRepeatsBedFile) : file('NO_FILE')
    rawVcf = runFreebayes(gatkResults.bamTuple, reorderFastaResults, genomeRepeatsBed)
    freebayesResults = filterAndSplitVcf(rawVcf)

    makeIndelTSV(freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi, consensusVcfGz, consensusVcfGzTbi ->
        tuple(sampleName, indelsVcfGz)
    })

    coverageBedResults = makeCoverageBed(gatkResults.bamTuple)

    makeConsensusFromVcfAndBed(
        freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi, consensusVcfGz, consensusVcfGzTbi ->
            tuple(sampleName, consensusVcfGz, consensusVcfGzTbi)
        }.join(coverageBedResults),
        reorderFastaResults
    )


    genomecovResults = genomecov(gatkResults.bamTuple, reorderFastaResults)

    bedgraphToBigWigResults = bedGraphToBigWig(reorderFastaResults, genomecovResults)

    if (runCnv) {
      sortForCountingResults = sortForCounting(gatkResults.bamTuple)

      htseqCountResults = htseqCount(sortForCountingResults, params.gtfFile)

      calculateTPMResults = calculateTPM(htseqCountResults, params.footprintFile)

      calculatePloidyAndGeneCNV(calculateTPMResults, params.footprintFile, params.ploidy, params.geneSourceIdOrthologFile, params.chrsForCalcFile)
    }

    makeWindowFileResults = makeWindowFile(reorderFastaResults, params.winLen)

    if (runCnv) {
      bedtoolsWindowedResults =  bedtoolsWindowed(makeWindowFileResults, gatkResults.bamTuple)

      normaliseCoverageResults = normaliseCoverage(bedtoolsWindowedResults, params.chrsForCalcFile, params.ploidy)

      normaliseCoverageToBigWigResults = normaliseCoverageToBigWig(reorderFastaResults, normaliseCoverageResults)
    }

    // CONVERT bed to bw here

    
    makeSnpDensityResults = makeSnpDensity(
        freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi, consensusVcfGz, consensusVcfGzTbi ->
            tuple(sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi)
        },
        makeWindowFileResults
    )

    makeDensityBigwigsResults = makeDensityBigwigs(makeSnpDensityResults, reorderFastaResults)

    if (params.ploidy != 1) {

      snpsVcf = freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi, consensusVcfGz, consensusVcfGzTbi ->
        tuple(sampleName, snpsVcfGz, snpsVcfGzTbi)
      }

      convertedSnpsVcf = convertFreebayesToVarscanFormat(snpsVcf)

      getHeterozygousSNPsResults = getHeterozygousSNPs(convertedSnpsVcf)

      makeHeterozygousDensityBedResults = makeHeterozygousDensityBed(makeWindowFileResults, getHeterozygousSNPsResults)

      makeHeterozygousDensityBigwig(makeHeterozygousDensityBedResults, reorderFastaResults)
    }

    // Alignment statistics: run samtools, bedtools, and raw read count in parallel,
    // then merge all per-sample results into a single published TSV
    rawReadCountResults = rawReadCount(files_only_qch)
    samtoolsStatsResults = samtoolsStats(gatkResults.bamTuple)
    bedtoolsGenomecovStatsResults = bedtoolsGenomecovStats(gatkResults.bamTuple, reorderFastaResults)
    makeAlignmentStats(
      samtoolsStatsResults
        .join(bedtoolsGenomecovStatsResults)
        .join(rawReadCountResults)
    )

}
