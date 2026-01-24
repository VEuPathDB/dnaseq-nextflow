#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Preprocessing
include { downloadBAMFromEBI } from '../modules/preprocessing.nf'
include { downloadFiles } from '../modules/preprocessing.nf'
include { fastqc } from '../modules/preprocessing.nf'
include { fastqc_check } from '../modules/preprocessing.nf'
include { trimmomatic } from '../modules/preprocessing.nf'

// Alignment
include { hisat2Index } from '../modules/alignment.nf'
include { hisat2 } from '../modules/alignment.nf'
include { reorderFasta } from '../modules/alignment.nf'
include { subsample } from '../modules/alignment.nf'
include { picard } from '../modules/alignment.nf'
include { gatk } from '../modules/alignment.nf'

// SNP
include { mpileup } from '../modules/snp.nf'
include { varscan } from '../modules/snp.nf'
include { concatSnpsAndIndels } from '../modules/snp.nf'
include { makeCombinedVarscanIndex } from '../modules/snp.nf'
include { filterIndels } from '../modules/snp.nf'
include { makeIndelTSV } from '../modules/snp.nf'
include { mergeVcfs } from '../modules/snp.nf'
include { makeMergedVarscanIndex } from '../modules/snp.nf'
include { bcftoolsConsensus } from '../modules/snp.nf'
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

    hisat2IndexResults = hisat2Index(genome_fasta_file, params.fromBAM, params.createIndex)

    if(!params.local && !params.fromBAM) {

        downloadFilesResults = downloadFiles(samples_qch)

        fastqcResults = fastqc(downloadFilesResults.files, params.fromBAM)

        fastqc_checkResults = fastqc_check(downloadFilesResults.files.join(fastqcResults), params.fromBAM)

        trimmomaticResults = trimmomatic(downloadFilesResults.files.join(fastqc_checkResults), params.fromBAM, downloadFilesResults.isPaired)

        hisat2Results = hisat2(downloadFilesResults.files.join(fastqc_checkResults).join(trimmomaticResults), hisat2IndexResults.genome_index_name, hisat2IndexResults.ht2_files, params.fromBAM, downloadFilesResults.isPaired)
    }

    else if(!params.local && params.fromBAM) {

        files = downloadBAMFromEBI(samples_qch)

        fastqcResults = fastqc(files, params.fromBAM)

        fastqc_checkResults = fastqc_check(files.join(fastqcResults), params.fromBAM)

        trimmomaticResults = trimmomatic(files.join(fastqc_checkResults), params.fromBAM, 'NA')

        hisat2Results = hisat2(files.join(fastqc_checkResults).join(trimmomaticResults), hisat2IndexResults.genome_index_name, hisat2IndexResults.ht2_files, params.fromBAM, 'NA')

    }

    else {

        fastqcResults = fastqc(samples_qch, params.fromBAM)

        fastqc_checkResults = fastqc_check(samples_qch.join(fastqcResults), params.fromBAM)

        trimmomaticResults = trimmomatic(samples_qch.join(fastqc_checkResults), params.fromBAM, params.isPaired)

        hisat2Results = hisat2(samples_qch.join(fastqc_checkResults).join(trimmomaticResults), hisat2IndexResults.genome_index_name, hisat2IndexResults.ht2_files, params.fromBAM, params.isPaired)

    }

    reorderFastaResults = reorderFasta(hisat2Results.first(), genome_fasta_file)

    subsampleResults = subsample(hisat2Results)

    picardResults = picard(reorderFastaResults, subsampleResults)

    gatkResults = gatk(reorderFastaResults, picardResults.bam_and_dict )

    mpileupResults = mpileup(gatkResults, reorderFastaResults)

    varscanResults = varscan(gatkResults.join(mpileupResults), reorderFastaResults)

    concatSnpsAndIndelsResults = concatSnpsAndIndels(varscanResults.vcf_files)

    makeCombinedVarscanIndexResults = makeCombinedVarscanIndex(concatSnpsAndIndelsResults)

    filterIndelsResults = filterIndels(makeCombinedVarscanIndexResults)

    makeIndelTSV(filterIndelsResults)

    // NOTE:  Must ensure the order here is consistent for the vcf files and their indexes;  the lists of paths are each sorted
    mergeVcfsResults = mergeVcfs(makeCombinedVarscanIndexResults.count(), makeCombinedVarscanIndexResults.map{ tuple it[1], it[2], "key" }.groupTuple(by: 2, sort: { a, b -> a <=> b } ))

    makeMergedVarscanIndexResults = makeMergedVarscanIndex(mergeVcfsResults)

    bcftoolsConsensusResults = bcftoolsConsensus(makeCombinedVarscanIndexResults, reorderFastaResults)

    addSampleToDefline(bcftoolsConsensusResults)

    genomecovResults = genomecov(gatkResults, reorderFastaResults)

    bedgraphToBigWigResults = bedGraphToBigWig(reorderFastaResults, genomecovResults)

    sortForCountingResults = sortForCounting(gatkResults)

    htseqCountResults = htseqCount(sortForCountingResults, params.gtfFile)

    calculateTPMResults = calculateTPM(htseqCountResults, params.footprintFile)

    calculatePloidyAndGeneCNV(calculateTPMResults, params.footprintFile, params.ploidy, params.taxonId, params.geneSourceIdOrthologFile, params.chrsForCalcFile)

    makeWindowFileResults = makeWindowFile(reorderFastaResults, params.winLen)

    bedtoolsWindowedResults =  bedtoolsWindowed(makeWindowFileResults, gatkResults)

    normaliseCoverageResults = normaliseCoverage(bedtoolsWindowedResults.join(picardResults.metrics))

    makeSnpDensityResults = makeSnpDensity(varscanResults.vcf_files, makeWindowFileResults)

    makeDensityBigwigsResults = makeDensityBigwigs(makeSnpDensityResults, reorderFastaResults)

    if (params.ploidy != 1) {

      getHeterozygousSNPsResults = getHeterozygousSNPs(varscanResults.vcf_files)

      makeHeterozygousDensityBedResults = makeHeterozygousDensityBed(makeWindowFileResults, getHeterozygousSNPsResults)

      makeHeterozygousDensityBigwig(makeHeterozygousDensityBedResults, reorderFastaResults)
    }

}
