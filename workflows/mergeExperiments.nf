#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { mergeVcfs }          from '../modules/mergeExperiments.nf'
include { mergeCoverageBeds }  from '../modules/mergeExperiments.nf'
include { makeGenomicIndelDb } from '../modules/mergeExperiments.nf'
include { makeCodingData }     from '../modules/mergeExperiments.nf'
include { processSeqVars }     from '../modules/mergeExperiments.nf'
include { snpEff }             from '../modules/mergeExperiments.nf'

workflow me {

  take:
    fastas_qch
    vcfs_qch
    indels_qch
    coverages_qch

  main:

    combinedIndels = indels_qch.collectFile(name: 'indel.tsv')

    genomicIndelDb = makeGenomicIndelDb(combinedIndels)

    allFastas = fastas_qch.collect()

    allVcfs       = vcfs_qch.collect()
    allVcfsBranch = allVcfs.branch { single: it.size() == 1; multiple: true }
    mergedVcf     = allVcfsBranch.single.map { it[0] }
                      .mix(mergeVcfs(allVcfsBranch.multiple))

    coverageTsv = mergeCoverageBeds(coverages_qch.collect())

    codingData = makeCodingData(allFastas, genomicIndelDb, params.gtfFile, params.genomeFastaFile)

    processSeqVarsResults = processSeqVars(
      mergedVcf,
      params.vcfCacheFile,
      params.undoneStrains,
      params.reference_strain,
      codingData.codingSequencesDb,
      codingData.codingIndelsDb,
      params.gtfFile,
      coverageTsv
    )

    snpEff(processSeqVarsResults.outputVcf, params.gtfFile, params.genomeFastaFile)

}
