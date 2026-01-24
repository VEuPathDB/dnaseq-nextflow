#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { writeIndelConfigFile } from '../modules/loadSingleExperiment.nf'
include { writePloidyConfigFile } from '../modules/loadSingleExperiment.nf'
include { writeCNVConfigFile } from '../modules/loadSingleExperiment.nf'
include { loadIndels } from '../modules/loadSingleExperiment.nf'
include { loadPloidy } from '../modules/loadSingleExperiment.nf'
include { loadGeneCNV } from '../modules/loadSingleExperiment.nf'

workflow ls {

  take:
    indels_qch
    ploidy_qch
    cnv_qch

  main:
    writeIndelConfigFileResults = writeIndelConfigFile(indels_qch)
    writePloidyConfigFileResults = writePloidyConfigFile(ploidy_qch)
    writeCNVConfigFileResults = writeCNVConfigFile(cnv_qch)
    loadIndels(writeIndelConfigFileResults, params.extDbRlsSpec)
    loadPloidy(writePloidyConfigFileResults, params.extDbRlsSpec)
    loadGeneCNV(writeCNVConfigFileResults, params.extDbRlsSpec)
}
