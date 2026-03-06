#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { mergeVcfs } from '../modules/mergeExperiments.nf'
include { makeGenomicIndelDb } from '../modules/mergeExperiments.nf'
include { makeCodingData } from '../modules/mergeExperiments.nf'
include { processSeqVars } from '../modules/mergeExperiments.nf'
include { snpEff } from '../modules/mergeExperiments.nf'

workflow me {

  take:
    fastas_qch
    vcfs_qch
    gvcfs_qch
    indels_qch

  main:

    combinedIndels = indels_qch.collectFile(name: 'indel.tsv')
    allFastas      = fastas_qch.collect()

    genomicIndelDb = makeGenomicIndelDb(combinedIndels)

    codingData = makeCodingData(allFastas, genomicIndelDb, params.gtfFile, params.genomeFastaFile)

    allvcfs = vcfs_qch.collect()

    mergedVcf = mergeVcfs(allvcfs)

    processSeqVarsResults = processSeqVars(mergedVcf, params.vcfCacheFile, params.undoneStrains, params.reference_strain, codingData.codingSequencesDb, codingData.codingIndelsDb, params.gtfFile)

    snpEff(mergedVcf, params.gtfFile, params.genomeFastaFile)

}
