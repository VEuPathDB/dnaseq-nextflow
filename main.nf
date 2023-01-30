#!/usr/bin/env nextflow
import nextflow.splitter.CsvSplitter
nextflow.enable.dsl=2

def fetchRunAccessions( tsv ) {
    def splitter = new CsvSplitter().options( header:true, sep:'\t' )
    def reader = new BufferedReader( new FileReader( tsv ) )
    splitter.parseHeader( reader )
    List<String> run_accessions = []
    Map<String,String> row
    while( row = splitter.fetchRecord( reader ) ) {
       run_accessions.add( row['run_accession'] )
    }
    return run_accessions
}

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
    input = fetchRunAccessions(params.input)
    samples_qch = Channel.fromList(input)
  }

  else {
    samples_qch = Channel.fromPath([params.input + '/**/*.fastq', params.input + '/**/*.fastq.gz', params.input + '/**/*.fq.gz'])
                                                                                             .map { file -> tuple(file.baseName, [file]) }
  }

  ps(samples_qch)

}


//---------------------------------------------------------------
// mergeExperiments
//---------------------------------------------------------------

workflow mergeExperiments {

  if(params.inputDir) {
    fastas_qch = Channel.fromPath(params.inputDir + '/*.fa.gz')
    vcfs_qch = Channel.fromPath(params.inputDir + '/result.vcf.gz')
    indels_qch = Channel.fromPath(params.inputDir + '/*.indel.tsv')
    coverage_qch = Channel.fromPath(params.varscanFilePath + '/*.coverage.txt')
  }

  else {
    throw new Exception("Missing parameter params.inputDir")
  }
   
  me(fastas_qch, vcfs_qch, indels_qch, coverage_qch)

}

//---------------------------------------------------------------
// loadSingleExperiment
//---------------------------------------------------------------

workflow loadSingleExperiment {

  if(!params.input) {
    throw new Exception("Missing parameter params.input")
  }

  else {
    bam_qch = Channel.fromPath([params.input + '/*.bam'], checkIfExists: true)
    bw_qch = Channel.fromPath([params.input + '/*.bw'], checkIfExists: true)
  }

  ls(bam_qch, bw_qch)
  
}

//---------------------------------------------------------------
// loadCNV
//---------------------------------------------------------------

workflow loadCNV {

  if(!params.input) {
    throw new Exception("Missing parameter params.input")
  }

  else {
    tpm_qch = Channel.fromPath([params.input + '/CNVs/*.tpm'], checkIfExists: true).map { file -> tuple(file.baseName, [file]) }
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

