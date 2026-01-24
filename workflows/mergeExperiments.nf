#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { checkUniqueIds } from '../modules/mergeExperiments.nf'
include { mergeVcfs } from '../modules/mergeExperiments.nf'
include { makeSnpFile } from '../modules/mergeExperiments.nf'
include { processSeqVars } from '../modules/mergeExperiments.nf'
include { addFeatureIdsToVariation } from '../modules/mergeExperiments.nf'
include { insertVariation } from '../modules/mergeExperiments.nf'
include { insertProduct } from '../modules/mergeExperiments.nf'
include { insertAllele } from '../modules/mergeExperiments.nf'
include { snpEff } from '../modules/mergeExperiments.nf'

workflow me {

  take:

    fastas_qch
    vcfs_qch
    indels_qch
    coverage_qch
    bw_qch
    bam_qch

  main:

    bigwigs = bw_qch.collectFile(storeDir: params.webServicesDir)
    bams = bam_qch.collectFile(storeDir: params.webServicesDir)

    coverages = coverage_qch.collectFile(storeDir: params.varscan_directory)

    combinedFastagz = fastas_qch.collectFile(name: 'CombinedFasta.fa.gz')
    combinedIndels = indels_qch.collectFile(name: 'indel.tsv')

    checkUniqueIds(combinedFastagz)

    allvcfs = vcfs_qch.collect()

    mergedVcf = mergeVcfs(allvcfs)

    makeSnpFileResults = makeSnpFile(mergedVcf)

    processSeqVarsResults = processSeqVars(makeSnpFileResults.snpFile, params.cacheFile, params.undoneStrains, params.organism_abbrev, params.reference_strain, params.varscan_directory, params.genomeFastaFile, combinedFastagz, combinedIndels, params.gtfFile, coverages, bigwigs, bams)

    addFeatureIdsToVariationResults = addFeatureIdsToVariation(processSeqVarsResults.variationFile, params.gusConfig)

    insertVariation(params.extDbRlsSpec, addFeatureIdsToVariationResults)
    insertProduct(processSeqVarsResults.productFile)
    insertAllele(processSeqVarsResults.alleleFile)

    snpEff(mergedVcf, params.gtfFile, params.genomeFastaFile)

}
