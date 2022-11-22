#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//---------------------------------------------------------------
// Including Workflows
//---------------------------------------------------------------

include { ps } from './modules/processSingleExperiment.nf'
include { ls } from './modules/loadSingleExperiment.nf'
include { lc } from './modules/loadCNV.nf'
include { me } from './modules/mergeExperiments.nf'
include { tests } from './modules/runTests.nf'

//---------------------------------------------------------------
// processSingleExperiment
//---------------------------------------------------------------

workflow processSingleExperiment {

  if(!params.inputDir) {
    throw new Exception("Missing parameter params.inputDir")
  }
   
  if(!params.genomeFastaFile) {
    throw new Exception("Missing parameter params.genomeFastaFile")
  }
  
  if(!params.gtfFile) {
    throw new Exception("Missing parameter params.gtfFile")
  }
  
  if(!params.geneFootprintFile) {
    throw new Exception("Missing parameter params.geneFootprintFile")
  }
  
  if(!params.trimmomaticAdaptorsFile) {
    throw new Exception("Missing parameter params.trimmomaticAdaptorsFile")
  }
  
  if(params.fromBAM) {
    samples_qch = Channel.fromPath([params.inputDir + '/**/*.bam'])
                    .map { file -> tuple(file.baseName, [file]) }
  }
  
  else if(params.isPaired) {
    samples_qch = Channel.fromFilePairs([params.inputDir + '/**/*_{1,2}.fastq', params.inputDir + '/**/*_{1,2}.fastq.gz', params.inputDir + '/**/*_{1,2}.fq.gz'])
  }
  
  else {
    samples_qch = Channel.fromPath([params.inputDir + '/**/*.fastq', params.inputDir + '/**/*.fastq.gz', params.inputDir + '/**/*.fq.gz'])
                                                                                             .map { file -> tuple(file.baseName, [file]) }
  }

  ps(samples_qch)

}

//---------------------------------------------------------------
// mergeExperiments
//---------------------------------------------------------------

workflow mergeExperiments {

  if(params.fastaDir) {
    fastas_qch = Channel.fromPath(params.fastaDir + '*.fa.gz')
  }
  
  else {
    throw new Exception("Missing parameter params.fastaDir")
  }
   
  if(params.vcfDir) {
    vcfs_qch = Channel.fromPath(params.vcfDir + '*.vcf.gz')
  }
   
  else {
    throw new Exception("Missing parameter params.vcfDir")
  }

  me(fastas_qch, vcfs_qch)

}

//---------------------------------------------------------------
// loadSingleExperiment
//---------------------------------------------------------------

workflow loadSingleExperiment {

  if(!params.inputDir) {
    throw new Exception("Missing parameter params.inputDir")
  }

  else {
    indels_qch = Channel.fromPath([params.inputDir + '/*.indel.tsv'], checkIfExists: true)
    bam_qch = Channel.fromPath([params.inputDir + '/*.bam'], checkIfExists: true)
    bw_qch = Channel.fromPath([params.inputDir + '/*.bw'], checkIfExists: true)
  }

  ls(indels_qch, bam_qch, bw_qch)
  
}

//---------------------------------------------------------------
// loadCNV
//---------------------------------------------------------------

workflow loadCNV {

  if(!params.inputDir) {
    throw new Exception("Missing parameter params.inputDir")
  }

  else {
    tpm_qch = Channel.fromPath([params.inputDir + '/CNVs/*.tpm'], checkIfExists: true).map { file -> tuple(file.baseName, [file]) }
  }

  lc(tpm_qch)
  
}

//---------------------------------------------------------------
// runTests
//---------------------------------------------------------------

workflow runTests {

  if(params.testDir) {
    tests_qch = Channel.fromPath([params.testDir + '*.t'])
  }
  
  else {
    throw new Exception("Missing parameter params.testDir")
  }

  tests(tests_qch)
  
}

//---------------------------------------------------------------
// DEFAULT - processSingleExperiment
//---------------------------------------------------------------

workflow {

  if(!params.inputDir) {
    throw new Exception("Missing parameter params.inputDir")
  }
   
  if(!params.genomeFastaFile) {
    throw new Exception("Missing parameter params.genomeFastaFile")
  }
  
  if(!params.gtfFile) {
    throw new Exception("Missing parameter params.gtfFile")
  }
  
  if(!params.geneFootprintFile) {
    throw new Exception("Missing parameter params.geneFootprintFile")
  }
  
  if(!params.trimmomaticAdaptorsFile) {
    throw new Exception("Missing parameter params.trimmomaticAdaptorsFile")
  }
  
  if(params.fromBAM) {
    samples_qch = Channel.fromPath([params.inputDir + '/**/*.bam'])
                    .map { file -> tuple(file.baseName, [file]) }
  }
  
  else if(params.isPaired) {
    samples_qch = Channel.fromFilePairs([params.inputDir + '/**/*_{1,2}.fastq', params.inputDir + '/**/*_{1,2}.fastq.gz', params.inputDir + '/**/*_{1,2}.fq.gz'])
  }
  
  else {
    samples_qch = Channel.fromPath([params.inputDir + '/**/*.fastq', params.inputDir + '/**/*.fastq.gz', params.inputDir + '/**/*.fq.gz'])
                                                                                             .map { file -> tuple(file.baseName, [file]) }
  }

  ps(samples_qch)

}

