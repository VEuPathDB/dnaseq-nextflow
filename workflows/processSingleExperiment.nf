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

// SNP
include { freebayes } from '../modules/snp.nf'
include { concatSnpsAndIndels } from '../modules/snp.nf'
include { makeCombinedVariantIndex } from '../modules/snp.nf'
include { filterIndels } from '../modules/snp.nf'
include { makeIndelTSV } from '../modules/snp.nf'
include { mergeVcfs } from '../modules/snp.nf'
include { makeMergedVariantIndex } from '../modules/snp.nf'
include { bcftoolsConsensusAndMask } from '../modules/snp.nf'
include { addSampleToDefline } from '../modules/snp.nf'

// CNV
include { genomecov } from '../modules/cnv.nf'
include { bedGraphToBigWig } from '../modules/cnv.nf'
include { bedGraphToBigWig as normaliseCoverageToBigWig } from '../modules/cnv.nf'

include { sortForCounting } from '../modules/cnv.nf'
include { htseqCount } from '../modules/cnv.nf'
include { calculateTPM } from '../modules/cnv.nf'
include { makeWindowFile } from '../modules/cnv.nf'
include { bedtoolsWindowed } from '../modules/cnv.nf'
include { normaliseCoverage } from '../modules/cnv.nf'
include { makeSnpDensity } from '../modules/cnv.nf'
include { makeDensityBigwigs } from '../modules/cnv.nf'
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

    freebayesResults = freebayes(gatkResults.bamTuple, reorderFastaResults)

    concatSnpsAndIndelsResults = concatSnpsAndIndels(freebayesResults.vcf_files)

    makeCombinedVariantIndexResults = makeCombinedVariantIndex(concatSnpsAndIndelsResults)

    filterIndelsResults = filterIndels(makeCombinedVariantIndexResults)

    makeIndelTSV(filterIndelsResults)

    // NOTE:  Must ensure the order here is consistent for the vcf files and their indexes;  the lists of paths are each sorted
    mergeVcfsResults = mergeVcfs(makeCombinedVariantIndexResults.count(), makeCombinedVariantIndexResults.map{ tuple it[1], it[2], "key" }.groupTuple(by: 2, sort: { a, b -> a <=> b } ))

    makeMergedVariantIndexResults = makeMergedVariantIndex(mergeVcfsResults)

    bcftoolsConsensusAndMaskResults = bcftoolsConsensusAndMask(makeCombinedVariantIndexResults, reorderFastaResults, gatkResults.bamFiles.collect())

    addSampleToDefline(bcftoolsConsensusAndMaskResults)

    genomecovResults = genomecov(gatkResults.bamTuple, reorderFastaResults)

    bedgraphToBigWigResults = bedGraphToBigWig(reorderFastaResults, genomecovResults)

    sortForCountingResults = sortForCounting(gatkResults.bamTuple)

    htseqCountResults = htseqCount(sortForCountingResults, params.gtfFile)

    calculateTPMResults = calculateTPM(htseqCountResults, params.footprintFile)

    calculatePloidyAndGeneCNV(calculateTPMResults, params.footprintFile, params.ploidy, params.taxonId, params.geneSourceIdOrthologFile, params.chrsForCalcFile)

    makeWindowFileResults = makeWindowFile(reorderFastaResults, params.winLen)

    bedtoolsWindowedResults =  bedtoolsWindowed(makeWindowFileResults, gatkResults.bamTuple)

    normaliseCoverageResults = normaliseCoverage(bedtoolsWindowedResults.join(picardResults.metrics))

    normaliseCoverageToBigWigResults = normaliseCoverageToBigWig(reorderFastaResults, normaliseCoverageResults)

    // CONVERT bed to bw here

    
    makeSnpDensityResults = makeSnpDensity(freebayesResults.vcf_files, makeWindowFileResults)

    makeDensityBigwigsResults = makeDensityBigwigs(makeSnpDensityResults, reorderFastaResults)

    if (params.ploidy != 1) {

      getHeterozygousSNPsResults = getHeterozygousSNPs(freebayesResults.vcf_files)

      makeHeterozygousDensityBedResults = makeHeterozygousDensityBed(makeWindowFileResults, getHeterozygousSNPsResults)

      makeHeterozygousDensityBigwig(makeHeterozygousDensityBedResults, reorderFastaResults)
    }

}
