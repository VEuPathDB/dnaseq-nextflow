#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//---------------------------------------------------------------
// Including Workflows
//---------------------------------------------------------------

include { ps } from './modules/processSingleExperiment.nf'
include { me } from './modules/mergeExperiments.nf'
include { ls } from './modules/loadSingleExperiment.nf'
include { wf } from './modules/wholeWorkflow.nf'
include { tests } from './modules/runTests.nf'

//---------------------------------------------------------------
// processSingleExperiment
//---------------------------------------------------------------

workflow processSingleExperiment {

  if(!params.input) {
    throw new Exception("Missing parameter params.input")
  }
   
  if(!params.genomeFastaFile) {
    throw new Exception("Missing parameter params.genomeFastaFile")
  }
  
  if(!params.gtfFile) {
    throw new Exception("Missing parameter params.gtfFile")
  }
  
  if(!params.footprintFile) {
    throw new Exception("Missing parameter params.footprintFile")
  }
  
  if(!params.trimmomaticAdaptorsFile) {
    throw new Exception("Missing parameter params.trimmomaticAdaptorsFile")
  }
  
  if(params.fromBAM && params.local) {
    samples_qch = Channel.fromPath([params.input + '/**/*.bam'])
                    .map { file -> tuple(file.baseName, [file]) }
  }

  else if(params.isPaired && params.local) {
    samples_qch = Channel.fromFilePairs([params.input + '/**/*_{1,2}.fastq', params.input + '/**/*_{1,2}.fastq.gz', params.input + '/**/*_{1,2}.fq.gz'])
  }

  else if(!params.local) {
    samples_qch = Channel.fromPath(params.input).splitCsv(sep: '\t')
  }

  else {
    samples_qch = Channel.fromPath([params.input + '/**/*.fastq', params.input + '/**/*.fastq.gz', params.input + '/**/*.fq.gz'])
                                                                                             .map { file -> tuple(file.baseName, [file]) }
  }

  ps(samples_qch)

}

//---------------------------------------------------------------
// loadSingleExperiment
//---------------------------------------------------------------

workflow loadSingleExperiment {

  if(params.indelDir) {
    indels_qch = Channel.fromPath(params.indelDir + '/*.indel.tsv')
  }

  else {
    throw new Exception("Missing parameter params.indelDir")
  }
   
  ls(indels_qch)
}

//---------------------------------------------------------------
// mergeExperiments
//---------------------------------------------------------------

workflow mergeExperiments {

  if(params.inputDir) {
    fastas_qch = Channel.fromPath(params.inputDir + '/*.fa.gz')
    vcfs_qch = Channel.fromPath(params.inputDir + '/result.vcf.gz')
    indels_qch = Channel.fromPath(params.inputDir + '/*.indel.tsv')
    bam_qch = Channel.fromPath(params.inputDir + '/*.bam')
    bw_qch = Channel.fromPath(params.inputDir + '/*.bw')
    coverage_qch = Channel.fromPath(params.varscanFilePath + '/*.coverage.txt')
  }

  else {
    throw new Exception("Missing parameter params.inputDir")
  }
   
  me(fastas_qch, vcfs_qch, indels_qch, coverage_qch, bam_qch, bw_qch)

}


//---------------------------------------------------------------
// wholeWorkflow
//---------------------------------------------------------------

workflow wholeWorkflow {

  if(!params.input) {
    throw new Exception("Missing parameter params.input")
  }
   
  if(!params.genomeFastaFile) {
    throw new Exception("Missing parameter params.genomeFastaFile")
  }
  
  if(!params.gtfFile) {
    throw new Exception("Missing parameter params.gtfFile")
  }
  
  if(!params.footprintFile) {
    throw new Exception("Missing parameter params.footprintFile")
  }
  
  if(!params.trimmomaticAdaptorsFile) {
    throw new Exception("Missing parameter params.trimmomaticAdaptorsFile")
  }
  
  if(params.fromBAM && params.local) {
    samples_qch = Channel.fromPath([params.input + '/**/*.bam'])
                    .map { file -> tuple(file.baseName, [file]) }
  }

  else if(params.isPaired && params.local) {
    samples_qch = Channel.fromFilePairs([params.input + '/**/*_{1,2}.fastq', params.input + '/**/*_{1,2}.fastq.gz', params.input + '/**/*_{1,2}.fq.gz'])
  }

  else if(!params.local) {
    samples_qch = Channel.fromPath(params.input).splitCsv(sep: '\t')
  }

  else {
    samples_qch = Channel.fromPath([params.input + '/**/*.fastq', params.input + '/**/*.fastq.gz', params.input + '/**/*.fq.gz'])
                                                                                             .map { file -> tuple(file.baseName, [file]) }
  }

  ps(samples_qch)

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

  if(!params.input) {
    throw new Exception("Missing parameter params.input")
  }
   
  if(!params.genomeFastaFile) {
    throw new Exception("Missing parameter params.genomeFastaFile")
  }
  
  if(!params.gtfFile) {
    throw new Exception("Missing parameter params.gtfFile")
  }
  
  if(!params.footprintFile) {
    throw new Exception("Missing parameter params.footprintFile")
  }
  
  if(!params.trimmomaticAdaptorsFile) {
    throw new Exception("Missing parameter params.trimmomaticAdaptorsFile")
  }
  
  if(params.fromBAM) {
    samples_qch = Channel.fromPath([params.input + '/**/*.bam'])
                    .map { file -> tuple(file.baseName, [file]) }
  }
  
  else if(params.isPaired) {
    samples_qch = Channel.fromFilePairs([params.input + '/**/*_{1,2}.fastq', params.input + '/**/*_{1,2}.fastq.gz', params.input + '/**/*_{1,2}.fq.gz'])
  }
  
  else {
    samples_qch = Channel.fromPath([params.input + '/**/*.fastq', params.input + '/**/*.fastq.gz', params.input + '/**/*.fq.gz'])
                                                                                             .map { file -> tuple(file.baseName, [file]) }
  }

  ps(samples_qch)

}

