#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//---------------------------------------------------------------
// Including Workflows
//---------------------------------------------------------------

include { ps } from './workflows/processSingleExperiment.nf'
include { me } from './workflows/mergeExperiments.nf'
include { ls } from './workflows/loadSingleExperiment.nf'
include { tests } from './modules/runTests.nf'

//---------------------------------------------------------------
// processSingleExperiment
//---------------------------------------------------------------

workflow processSingleExperiment {

  if(!params.samplesheet) {
    throw new Exception("Missing parameter params.samplesheet")
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

  // Parse samplesheet and create input channel
  Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
      def sample_id = row.sample
      def fastq_1 = file(row.fastq_1, checkIfExists: true)
      def fastq_2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null
      def is_paired = row.fastq_2 ? true : false

      // Create file list based on pairing
      def files = fastq_2 ? [fastq_1, fastq_2] : [fastq_1]

      tuple(sample_id, files, is_paired)
    }
    .set { samples_qch }

  ps(samples_qch)

}


// TODO:  Need to sort out how to load into Postgres
//---------------------------------------------------------------
// loadSingleExperiment
//---------------------------------------------------------------

// workflow loadSingleExperiment {

//   if(params.inputDir) {
//     indels_qch = Channel.fromPath(params.inputDir + '/*.indel.tsv').map { file -> tuple(file.baseName, [file]) }
//     ploidy_qch = Channel.fromPath(params.inputDir + '/*_Ploidy.txt').map { file -> tuple(file.baseName, [file]) }
//     cnv_qch = Channel.fromPath(params.inputDir + '/_geneCNVs.txt').map { file -> tuple(file.baseName, [file]) }
//   }

//   else {
//     throw new Exception("Missing parameter params.indelDir")
//   }
   
//   ls(indels_qch,ploidy_qch,cnv_qch)
// }

//---------------------------------------------------------------
// mergeExperiments
//---------------------------------------------------------------

workflow mergeExperiments {
  fastas_qch    = Channel.fromPath(params.relativeConsensusFilePattern)
  vcfs_qch      = Channel.fromPath(params.vcfFiles)
  indels_qch    = Channel.fromPath(params.indelsFiles)
  coverages_qch = Channel.fromPath(params.coverageFiles)

  me(fastas_qch, vcfs_qch, indels_qch, coverages_qch)
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

  if(!params.samplesheet) {
    throw new Exception("Missing parameter params.samplesheet")
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

  // Parse samplesheet and create input channel
  Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
      def sample_id = row.sample
      def fastq_1 = file(row.fastq_1, checkIfExists: true)
      def fastq_2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null
      def is_paired = row.fastq_2 ? true : false

      // Create file list based on pairing
      def files = fastq_2 ? [fastq_1, fastq_2] : [fastq_1]

      tuple(sample_id, files, is_paired)
    }
    .set { samples_qch }

  ps(samples_qch)

}

