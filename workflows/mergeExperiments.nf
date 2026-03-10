#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { mergeVcfs; mergeVcfs as mergeGvcfs } from '../modules/mergeExperiments.nf'
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

    genomicIndelDb = makeGenomicIndelDb(combinedIndels)

    allFastas      = fastas_qch.collect()

    allvcfs  = vcfs_qch.collect().branch { single: it.size() == 1; multiple: true }
    allgvcfs = gvcfs_qch.collect().branch { single: it.size() == 1; multiple: true }

    mergedVcf  = allvcfs.single.map { it[0] }.mix(mergeVcfs(allvcfs.multiple))
    mergedGvcf = allgvcfs.single.map { it[0] }.mix(mergeGvcfs(allgvcfs.multiple))

    codingData = makeCodingData(allFastas, genomicIndelDb, params.gtfFile, params.genomeFastaFile)

    // processSeqVarsResults = processSeqVars(mergedVcf, params.vcfCacheFile, params.undoneStrains, params.reference_strain, codingData.codingSequencesDb, codingData.codingIndelsDb, params.gtfFile)

    // snpEff(mergedVcf, params.gtfFile, params.genomeFastaFile)

}
