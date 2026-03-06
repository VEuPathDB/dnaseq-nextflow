#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { checkUniqueIds } from '../modules/mergeExperiments.nf'
include { mergeVcfs } from '../modules/mergeExperiments.nf'
include { makeSnpFile } from '../modules/mergeExperiments.nf'
include { processSeqVars } from '../modules/mergeExperiments.nf'
include { snpEff } from '../modules/mergeExperiments.nf'

workflow me {

  take:
    fastas_qch
    vcfs_qch
    gvcfs_qch
    indels_qch

  main:

    combinedFastagz = fastas_qch.collectFile(name: 'CombinedFasta.fa.gz')
    combinedIndels = indels_qch.collectFile(name: 'indel.tsv')

    checkUniqueIds(combinedFastagz)

    allvcfs = vcfs_qch.collect()

    mergedVcf = mergeVcfs(allvcfs)

    makeSnpFileResults = makeSnpFile(mergedVcf)

    // Placeholder channels — wired to upstream transcript-prep process once it exists
    transcriptDb = Channel.of(params.transcriptDb)
    indelDb = Channel.of(params.indelDb)

    processSeqVarsResults = processSeqVars(makeSnpFileResults.snpFile, params.vcfCacheFile, params.undoneStrains, params.reference_strain, transcriptDb, indelDb, params.gtfFile)

    snpEff(mergedVcf, params.gtfFile, params.genomeFastaFile)

}
