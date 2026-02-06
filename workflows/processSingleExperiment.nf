#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Preprocessing
include { downloadBAMFromEBI } from '../modules/preprocessing.nf'
include { downloadFiles } from '../modules/preprocessing.nf'
include { fastqc } from '../modules/preprocessing.nf'
include { fastqc_check } from '../modules/preprocessing.nf'
include { trimmomatic } from '../modules/preprocessing.nf'

// Alignment
include { bwaIndex } from '../modules/alignment.nf'
include { bwaMem } from '../modules/alignment.nf'
include { reorderFasta } from '../modules/alignment.nf'
include { subsample } from '../modules/alignment.nf'
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

    bwaIndexResults = bwaIndex(genome_fasta_file, params.fromBAM, params.createIndex)

    if(!params.local && !params.fromBAM) {

        downloadFilesResults = downloadFiles(samples_qch)

        fastqcResults = fastqc(downloadFilesResults.files, params.fromBAM)

        fastqc_checkResults = fastqc_check(downloadFilesResults.files.join(fastqcResults), params.fromBAM)

        trimmomaticResults = trimmomatic(downloadFilesResults.files.join(fastqc_checkResults), params.fromBAM, downloadFilesResults.isPaired)

        bwaMemResults = bwaMem(downloadFilesResults.files.join(fastqc_checkResults).join(trimmomaticResults), bwaIndexResults.index_files, bwaIndexResults.genome_fasta, params.fromBAM, downloadFilesResults.isPaired)
    }

    else if(!params.local && params.fromBAM) {

        files = downloadBAMFromEBI(samples_qch)

        fastqcResults = fastqc(files, params.fromBAM)

        fastqc_checkResults = fastqc_check(files.join(fastqcResults), params.fromBAM)

        trimmomaticResults = trimmomatic(files.join(fastqc_checkResults), params.fromBAM, 'NA')

        bwaMemResults = bwaMem(files.join(fastqc_checkResults).join(trimmomaticResults), bwaIndexResults.index_files, bwaIndexResults.genome_fasta, params.fromBAM, 'NA')

    }

    else {

        fastqcResults = fastqc(samples_qch, params.fromBAM)

        fastqc_checkResults = fastqc_check(samples_qch.join(fastqcResults), params.fromBAM)

        trimmomaticResults = trimmomatic(samples_qch.join(fastqc_checkResults), params.fromBAM, params.isPaired)

        bwaMemResults = bwaMem(samples_qch.join(fastqc_checkResults).join(trimmomaticResults), bwaIndexResults.index_files, bwaIndexResults.genome_fasta, params.fromBAM, params.isPaired)

    }

    reorderFastaResults = reorderFasta(bwaMemResults.first(), genome_fasta_file)

    subsampleResults = subsample(bwaMemResults)

    picardResults = picard(reorderFastaResults, subsampleResults)

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

    makeSnpDensityResults = makeSnpDensity(freebayesResults.vcf_files, makeWindowFileResults)

    makeDensityBigwigsResults = makeDensityBigwigs(makeSnpDensityResults, reorderFastaResults)

    if (params.ploidy != 1) {

      getHeterozygousSNPsResults = getHeterozygousSNPs(freebayesResults.vcf_files)

      makeHeterozygousDensityBedResults = makeHeterozygousDensityBed(makeWindowFileResults, getHeterozygousSNPsResults)

      makeHeterozygousDensityBigwig(makeHeterozygousDensityBedResults, reorderFastaResults)
    }

}
