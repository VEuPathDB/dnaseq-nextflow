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
include { samtoolsStats } from '../modules/alignment.nf'
include { bedtoolsGenomecovStats } from '../modules/alignment.nf'
include { mergeAlignmentStats } from '../modules/alignment.nf'

// SNP
include { runFreebayes } from '../modules/snp.nf'
include { filterAndSplitVcf } from '../modules/snp.nf'
include { makeIndelTSV } from '../modules/snp.nf'
include { sanitizeVcf } from '../modules/snp.nf'
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

    rawVcf = runFreebayes(gatkResults.bamTuple, reorderFastaResults)
    freebayesResults = filterAndSplitVcf(rawVcf)

    filteredVcf = freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, vcfGz, vcfGzTbi)
    }

    // Feed the indels VCF produced by freebayes directly, bypassing the former filterIndels step
    makeIndelTSV(freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, indelsVcfGz)
    }).collectFile(name: 'indels.tsv', storeDir: params.outputDir)

    sanitizeVcfResults = sanitizeVcf(filteredVcf)

    coverageBedResults = makeCoverageBed(gatkResults.bamTuple)

    makeConsensusFromVcfAndBed(sanitizeVcfResults.join(coverageBedResults), reorderFastaResults)


    genomecovResults = genomecov(gatkResults.bamTuple, reorderFastaResults)

    bedgraphToBigWigResults = bedGraphToBigWig(reorderFastaResults, genomecovResults)

    sortForCountingResults = sortForCounting(gatkResults.bamTuple)

    htseqCountResults = htseqCount(sortForCountingResults, params.gtfFile)

    calculateTPMResults = calculateTPM(htseqCountResults, params.footprintFile)

    calculatePloidyAndGeneCNV(calculateTPMResults, params.footprintFile, params.ploidy, params.geneSourceIdOrthologFile, params.chrsForCalcFile)

    makeWindowFileResults = makeWindowFile(reorderFastaResults, params.winLen)

    bedtoolsWindowedResults =  bedtoolsWindowed(makeWindowFileResults, gatkResults.bamTuple)

    normaliseCoverageResults = normaliseCoverage(bedtoolsWindowedResults, params.chrsForCalcFile, params.ploidy)

    normaliseCoverageToBigWigResults = normaliseCoverageToBigWig(reorderFastaResults, normaliseCoverageResults)

    // CONVERT bed to bw here

    
    makeSnpDensityResults = makeSnpDensity(freebayesResults.vcf_files, makeWindowFileResults)

    makeDensityBigwigsResults = makeDensityBigwigs(makeSnpDensityResults, reorderFastaResults)

    if (params.ploidy != 1) {

      snpsVcf = freebayesResults.vcf_files.map { sampleName, vcfGz, vcfGzTbi, snpsVcfGz, snpsVcfGzTbi, indelsVcfGz, indelsVcfGzTbi ->
        tuple(sampleName, snpsVcfGz, snpsVcfGzTbi)
      }

      convertedSnpsVcf = convertFreebayesToVarscanFormat(snpsVcf)

      getHeterozygousSNPsResults = getHeterozygousSNPs(convertedSnpsVcf)

      makeHeterozygousDensityBedResults = makeHeterozygousDensityBed(makeWindowFileResults, getHeterozygousSNPsResults)

      makeHeterozygousDensityBigwig(makeHeterozygousDensityBedResults, reorderFastaResults)
    }

    // Alignment statistics: run samtools and bedtools in parallel on the final GATK BAM,
    // then merge all per-sample results into a single published TSV
    samtoolsStatsResults = samtoolsStats(gatkResults.bamTuple)
    bedtoolsGenomecovStatsResults = bedtoolsGenomecovStats(gatkResults.bamTuple, reorderFastaResults)
    mergeAlignmentStats(
      samtoolsStatsResults.map { sampleName, tsv -> tsv }.collect(),
      bedtoolsGenomecovStatsResults.map { sampleName, tsv -> tsv }.collect()
    )

}
