#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { mergeVcfs }          from '../modules/mergeExperiments.nf'
include { mergeCoverageBeds }  from '../modules/mergeExperiments.nf'
include { makeGenomicIndelDb } from '../modules/mergeExperiments.nf'
include { makeCodingData }     from '../modules/mergeExperiments.nf'
include { processSeqVars }     from '../modules/mergeExperiments.nf'
include { snpEff }                 from '../modules/mergeExperiments.nf'
include { parseSnpEffAnnotations } from '../modules/mergeExperiments.nf'

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

    // Arity is mergeVcfs' problem, not the DAG's: bin/mergeVcfs.sh normalizes
    // every input and merges only when there are 2+. Branching here previously
    // let the n=1 case skip normalization entirely.
    mergedVcf = mergeVcfs(vcfs_qch.collect())

    coverageTsv = mergeCoverageBeds(coverages_qch.collect())

    codingData = makeCodingData(allFastas, genomicIndelDb, params.gtfFile, params.genomeFastaFile)

    processSeqVarsResults = processSeqVars(
      mergedVcf,
      params.cacheFile,
      params.reference_strain,
      codingData.codingSequencesDb,
      codingData.codingIndelsDb,
      params.gtfFile,
      coverageTsv
    )

    snpEff(processSeqVarsResults.outputVcf, params.gtfFile, params.genomeFastaFile)

    parseSnpEffAnnotations(snpEff.out.annotatedVcf)

}
