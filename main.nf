#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//---------------------------------------------------------------
// Param Checking 
//---------------------------------------------------------------

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

//---------------------------------------------------------------
// Includes
//---------------------------------------------------------------

include { dnaseqAnalysis } from './modules/dnaseqAnalysis.nf'

//---------------------------------------------------------------
// Main Workflow
//---------------------------------------------------------------

workflow {

  dnaseqAnalysis(samples_qch)

}